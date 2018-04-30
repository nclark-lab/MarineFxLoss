#!/usr/bin/perl -w
use strict;
use PerlSubs;

# The script will stream through all genes in the UCSC knownCanonical.exonAA file and score stop codons.
# Stop codons are scored as potential pseudogene calls if they meet strict criteria.
# Also reports the proportion of considered gene sequence that is gap/missing for each species

### Pseudogene call criteria ####
my $firstExon_th = 0.10; #A. threshold proportion of total protein below which a first exon is completely excluded from analysis
my $firstExon_exclude = 20; #B. number of amino acids (codons) to exclude from the first exon if first exon is less then B% of total protein
my $lastExon_th = 0.10;
my $lastExon_exclude = 20;
my $maxPropExon = 0.25; # threshold proportion of exon that is made up of gaps above which to exclude all stops within that exon
#################################

my %stops; # stop codon tabulation by {gene} and {species}
my %proportion; # proportion of gap characters by {gene} and {species}
my ($g , $sp , $ex , $t_ex) = ('','','',''); # store gene name, species, exon # and total exons number
my %seq;
my $count = 1;
my $first_time = 1;
my $stopoutfile = 'LatestStopsExcludeExon25.tsv'; #output stops file
my $missoutfile = 'AAMissingExcludeExon25.tsv'; #output missingness file
my $errfile = 'LatestStopsExcludeExon25.err'; #error file

open (FILE , 'knownCanonical.exonAA.fa');
#open (FILE , 'testAA.fa');
#open (FILE , 'minitestAA.fa');
#open (FILE , 'ExonsAllSLC39s.fa');
#open (FILE , 'test.fa'); #input file

open(my $efh, '>', $errfile) or die "Could not open file '$errfile' $!";

while (<FILE>) {
	chomp;
	next if /^\s*$/; # skip pure whitespace lines

	# Break fasta header line into information components	
	if (/^>(\S+)_(\S+)_(\d+)_(\d+)\s/) {
		if ($first_time) {
			($g , $sp , $ex , $t_ex) = ($1 , $2 , $3 , $4);
			$first_time = 0;
		} else {
		
			unless ($g eq $1) {
				# score stop codons "Z"
				die "!!! Warning, exons don't add up. $g	$ex	$t_ex !!!\n" unless ($ex == $t_ex);
				print $efh "$count	Scoring stops in $g\n";
				
				#tabulate length of entire protein and first and last exons in human (hg19)
				$seq{'hg19'}->[0] = '';
				my $human_seq = join '', @{$seq{'hg19'}};
				my $human_len = length $human_seq;
				my $firstExon_prop = (length $seq{'hg19'}->[1]) / $human_len;
				my $lastExon_prop = (length $seq{'hg19'}->[$t_ex]) / $human_len;
				
				# step through species
				foreach my $this_sp (keys %seq) {
					# initiate stop codon and gap counters for gene and species
					$stops{$g}->{$this_sp} = 0;
					my $gene_gaps = 0;
                                        my $gene_len = 0;

					# step through exons (1 .. $t_ex)
					for (my $this_ex=1; $this_ex<=$t_ex; $this_ex++) {
					
						# score first exon (1)
						if ($this_ex == 1) {
							if ( $firstExon_prop < $firstExon_th ) { # skip exon 1 if it is less than threshold proportion of entire protein
								next;
							} else { # otherwise, first exon encodes large proportion of protein, so just ignore first $firstExon_exclude residues.
								#If this is the only exon...
								if ($t_ex == 1) {
									# next if length is <= sum of 2 first and last "excludes", next
									next if (length $seq{$this_sp}->[$this_ex]) <= ($firstExon_exclude + $lastExon_exclude);
									# trim front and back
									my $this_seq = substr $seq{$this_sp}->[$this_ex] , $firstExon_exclude, -$lastExon_exclude;
									my $n_Z = $this_seq =~ s/Z/X/g;
                                                                        my $seq_len_with_gaps = length $this_seq;
									my $totgaps = $this_seq =~ s/-//g;
									if ($totgaps / $seq_len_with_gaps < $maxPropExon) {
										$stops{$g}->{$this_sp} += $n_Z;
                                                                                $gene_gaps += $totgaps;
									} else {
                                                                                $gene_gaps += $seq_len_with_gaps;
                                                                        }
                                                                        $gene_len += $seq_len_with_gaps;
								} else {
									next if (length $seq{$this_sp}->[$this_ex]) <= $firstExon_exclude;
									my $this_seq = substr $seq{$this_sp}->[$this_ex] , $firstExon_exclude;
									my $n_Z = $this_seq =~ s/Z/X/g;
                                                                        my $seq_len_with_gaps = length $this_seq;
                                                                        my $totgaps = $this_seq =~ s/-//g;
                                                                        if ($totgaps / $seq_len_with_gaps < $maxPropExon) {
										$stops{$g}->{$this_sp} += $n_Z;
										$gene_gaps += $totgaps;
									} else {
										$gene_gaps += $seq_len_with_gaps;
									}
									$gene_len += $seq_len_with_gaps;
								}
							}
						
						# score last exon ($t_ex)
						} elsif ($this_ex == $t_ex) {
							if ( $lastExon_prop < $lastExon_th ) { # skip last exon if it is less than threshold proportion of entire protein
								next;
							} else { # otherwise, last exon encodes large proportion of protein, so just ignore first $lastExon_exclude residues.
								next if (length $seq{$this_sp}->[$this_ex]) <= $lastExon_exclude;
								my $this_seq = substr $seq{$this_sp}->[$this_ex] , 0 , -$lastExon_exclude;
								my $n_Z = $this_seq =~ s/Z/X/g;
                                                                my $seq_len_with_gaps = length $this_seq;
								my $totgaps = $this_seq =~ s/-//g;
								if ($totgaps / $seq_len_with_gaps  < $maxPropExon) {
									$stops{$g}->{$this_sp} += $n_Z;
                                                                        $gene_gaps += $totgaps;
								} else {
                                                                        $gene_gaps += $seq_len_with_gaps;
                                                                }
                                                                $gene_len += $seq_len_with_gaps;
							}
						
						# score internal exons (2 to $t_ex - 1) 
						} else {
							my $this_seq = $seq{$this_sp}->[$this_ex];
							my $n_Z = $this_seq =~ s/Z/X/g;
                                                        my $seq_len_with_gaps = length $this_seq;
							my $totgaps = $this_seq =~ s/-//g;
							if ($totgaps / $seq_len_with_gaps  < $maxPropExon) {
								$stops{$g}->{$this_sp} += $n_Z;
                                                                $gene_gaps += $totgaps;
                                                        } else {
                                                                $gene_gaps += $seq_len_with_gaps;
							}
                                                        $gene_len += $seq_len_with_gaps;
						}
					}

                                        # report proportion of gap characters for gene and species
                                        $proportion{$g}->{$this_sp} = $gene_gaps / $gene_len;				
                            }
				
				# Done with gene, reset seq and increment count.
				%seq = ();
				$count++;
			}

			($g , $sp , $ex , $t_ex) = ($1 , $2 , $3 , $4);
		}
	} else {
		$seq{$sp}->[$ex] = $_;	
	}
}
close FILE;
print $efh "Scans complete.\n";

## OUTPUT ##
# read in gene names
my %name;
foreach my $line (read_input_clean('names_nospaces.tsv')) { #removed spaces from names for post-processing
	my @e = split /\t/ , $line;
	$name{$e[0]} = $e[2];
}

# List all species
my @all_species = sort keys %seq;
open(my $sfh, '>', $stopoutfile) or die "Could not open file '$stopoutfile' $!";
print $sfh (join "\t" , 'gene', 'ucid', @all_species);
print $sfh "\n";

foreach my $gene_name (sort keys %stops) {
	if (defined $name{$gene_name}) {
		print $sfh "$name{$gene_name}\t$gene_name";
	} else { print "$gene_name\t$gene_name"; }
	foreach my $this_sp ( @all_species ) {
		print $sfh "\t$stops{$gene_name}->{$this_sp}";
	}
	print $sfh "\n";
}
close $sfh;

# Also write missingness proportion
open(my $mfh, '>', $missoutfile) or die "Could not open file '$missoutfile' $!";
print $mfh (join "\t" , 'gene', 'ucid', @all_species);
print $mfh "\n";

foreach my $gene_name (sort keys %proportion) {
        if (defined $name{$gene_name}) {
                print $mfh "$name{$gene_name}\t$gene_name";
        } else { print $mfh "$gene_name\t$gene_name"; }

        foreach my $this_sp ( @all_species ) {
                printf $mfh "\t%.4f" , $proportion{$gene_name}->{$this_sp};
        }
        print $mfh "\n";
}
close $mfh;
close $efh;

exit;

