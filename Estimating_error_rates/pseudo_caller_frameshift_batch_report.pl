#!/usr/bin/perl -w
use strict;
use PerlSubs;

# The script will stream through all genes in the UCSC knownCanonical.exonAA file and score stop codons.
# Stop codons are scored as potential pseudogene calls if they meet strict criteria.

### Pseudogene call criteria ####
my $firstExon_th = 0.10; #A. threshold proportion of total coding sequence below which a first exon is completely excluded from analysis
my $firstExon_exclude = 60; #B. number of nucleotides to exclude from the first exon if first exon is less then B% of total coding sequence
my $lastExon_th = 0.10;
my $lastExon_exclude = 60;
my $trim_junction = 3; # trim for all exonic sequence near splice site junctions
my $maxDel = 8; # largest permitted frameshifting deletion. Too large would confound with sequencing gaps.
my $proximity_th = 15; # proximity in basepairs between indepdent gap runs under which both are rejected.
my $maxPropExon = 0.25; # threshold proportion of exon that is made up of gaps above which to exclude all frameshifts within that exon
#################################

my %shifts; # frameshift tabulation by {gene} and {species}
my ($g , $sp , $ex , $t_ex) = ('','','',''); # store gene name, species, exon # and total exons number
my %seq;
my $count = 1;
my $first_time = 1;

# 58 mammals
my %species_include;
foreach (qw(hg19 panTro4 gorGor3 ponAbe2 nomLeu3 rheMac3 macFas5 papHam1 chlSab1 calJac3 saiBol1 otoGar3 tupChi1 speTri2 jacJac1 micOch1 criGri1 mesAur1 mm10 rn5 hetGla2 cavPor3 chiLan1 octDeg1 oryCun2 ochPri3 susScr3 vicPac2 camFer1 turTru2 orcOrc1 panHod1 bosTau7 oviAri3 capHir1 equCab2 cerSim1 felCat5 canFam3 musFur1 ailMel1 odoRosDiv1 lepWed1 pteAle1 pteVam1 myoDav1 myoLuc2 eptFus1 eriEur2 sorAra2 conCri1 loxAfr3 eleEdw1 triMan1 chrAsi1 echTel2 oryAfe1 dasNov3)) {
	$species_include{$_} = 1;
}

# read in gene names
my %name;
foreach my $line (read_input_clean('names_nospaces.tsv')) { #removed spaces from names for post-processing
	my @e = split /\t/ , $line;
	$name{$e[0]} = $e[2];
}

#open (FILE, 'ExonsAllSLC39s.fa');
#open (FILE , 'knownCanonical.exonNuc.fa');
#open (FILE , 'testNuc.fa');
open (FILE , 'knownCanonical.exonNuc.45genes.nosep.fa');
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
				# score frameshifts
				die "!!! Warning, exons don't add up. $g	$ex	$t_ex !!!\n" unless ($ex == $t_ex);
				print STDERR "$count	Scoring frameshifts in $g\n";
				
				#tabulate length of entire coding sequence and first and last exons in human (hg19)
				$seq{'hg19'}->[0] = '';
				my $human_seq = join '', @{$seq{'hg19'}};
				my $human_len = length $human_seq;
				my $firstExon_prop = (length $seq{'hg19'}->[1]) / $human_len;
				my $lastExon_prop = (length $seq{'hg19'}->[$t_ex]) / $human_len;
				
				# step through species
				foreach my $this_sp (sort keys %seq) {
					# initiate frameshift counter for gene and species
					$shifts{$g}->{$this_sp} = 0;
					my $L_report = 0;	#running length of species sequence up to previous exon
					
					# step through exons (1 .. $t_ex)
					for (my $this_ex=1; $this_ex<=$t_ex; $this_ex++) {
					
						# score first exon (1)
						if ($this_ex == 1) {
							if ( $firstExon_prop < $firstExon_th ) { # skip exon 1 if it is less than threshold proporiton of entire coding sequence
								next;
							} else { # otherwise, first exon encodes large proportion of coding sequence, so just ignore first $firstExon_exclude nucleotides.
								#If this is the only exon...
								if ($t_ex == 1) {
									# next if length is <= sum of 2 first and last "excludes", next
									next if (length $seq{$this_sp}->[$this_ex]) <= ($firstExon_exclude + $lastExon_exclude);
									# trim front and back
									my $this_seq = substr $seq{$this_sp}->[$this_ex] , $firstExon_exclude, -$lastExon_exclude;

									my ($n_Z , $starts , $stops) = score_frameshifts($this_seq);
									my @startsplus = map { $_ + $firstExon_exclude } @{$starts};
									my @stopsplus = map { $_ + $firstExon_exclude } @{$stops};
									if ($n_Z) { report($g , $this_sp , $seq{$this_sp}->[$this_ex] , $n_Z , \@startsplus , \@stopsplus); }
									$shifts{$g}->{$this_sp} += $n_Z;
								} else {
									next if (length $seq{$this_sp}->[$this_ex]) <= $firstExon_exclude;
									my $this_seq = substr $seq{$this_sp}->[$this_ex] , $firstExon_exclude , -$trim_junction;
									my ($n_Z , $starts , $stops) = score_frameshifts($this_seq);
									my @startsplus = map { $_ + $firstExon_exclude } @{$starts};
									my @stopsplus = map { $_ + $firstExon_exclude } @{$stops};
									if ($n_Z) { report($g , $this_sp , $seq{$this_sp}->[$this_ex] , $n_Z , \@startsplus , \@stopsplus); }
									$shifts{$g}->{$this_sp} += $n_Z;
								}
							}
						
							$L_report += length $seq{$this_sp}->[$this_ex];
							
						# score last exon ($t_ex)
						} elsif ($this_ex == $t_ex) {
							if ( $lastExon_prop < $lastExon_th ) { # skip last exon if it is less than threshold proporiton of entire protein
								next;
							} else { # otherwise, last exon encodes large proportion of protein, so just ignore first $lastExon_exclude residues.
								next if (length $seq{$this_sp}->[$this_ex]) <= $lastExon_exclude;
								my $this_seq = substr $seq{$this_sp}->[$this_ex] , $trim_junction , -$lastExon_exclude;
								my ($n_Z , $starts , $stops) = score_frameshifts($this_seq);
								my @startsplus = map { $_ + $trim_junction } @{$starts};
								my @stopsplus = map { $_ + $trim_junction } @{$stops};
								if ($n_Z) { report($g , $this_sp , $seq{$this_sp}->[$this_ex] , $n_Z , \@startsplus , \@stopsplus); }
								$shifts{$g}->{$this_sp} += $n_Z;
							}
						
							$L_report += length $seq{$this_sp}->[$this_ex];

						# score internal exons (2 to $t_ex - 1) 
						} else {
							my $this_seq = substr $seq{$this_sp}->[$this_ex] , $trim_junction , -$trim_junction ; # trim near splice sites
							my ($n_Z , $starts , $stops) = score_frameshifts($this_seq);
							my @startsplus = map { $_ + $trim_junction } @{$starts};
							my @stopsplus = map { $_ + $trim_junction } @{$stops};
							if ($n_Z) { report($g , $this_sp , $seq{$this_sp}->[$this_ex] , $n_Z , \@startsplus , \@stopsplus); }
							$shifts{$g}->{$this_sp} += $n_Z;
							
							$L_report += length $seq{$this_sp}->[$this_ex];

						}
					}
				}
				
				# Done with gene, reset seq and increment count.
				%seq = ();
				$count++;
			}

			($g , $sp , $ex , $t_ex) = ($1 , $2 , $3 , $4);
		}
	} else {
		if ($species_include{$sp}) {
			$seq{$sp}->[$ex] = $_;
		}
	}
}
close FILE;
print STDERR "Scans complete.\n";

## OUTPUT ##
# List all species
# my @all_species = sort keys %seq;
# print (join "\t" , 'gene', 'ucid', @all_species);
# print "\n";
# 
# foreach my $gene_name (sort keys %shifts) {
# 	if (defined $name{$gene_name}) {
# 		print "$name{$gene_name}\t$gene_name";
# 	} else { print "$gene_name\t$gene_name"; }
# 	foreach my $this_sp ( @all_species ) {
# 		print "\t$shifts{$gene_name}->{$this_sp}";
# 	}
# 	print "\n";
# }

exit;

sub report {
	my ($g , $this_sp , $this_seq , $n_Z , $starts , $stops) = @_;
	my $common_name = $g;
	if (defined $name{$g}) { $common_name = $name{$g}; }
	for (my $i=0; $i < $n_Z; $i++) {
		my $length = $stops->[$i] - $starts->[$i] + 1;
		my $this_start = $starts->[$i] - 1;
		my $this_stop = $stops->[$i] + 1;
		print STDOUT "$g	$common_name	$this_sp	frameshift	$length	$this_start	$this_stop	$this_seq\n";  ##############################
	}
}

sub score_frameshifts {
	# retrieve sequence
	my @seq = split '' , $_[0];
	
	# collect start and stop coordinates of gap character '-' runs
	my @start; my @stop;
	my $in = 0;
    	my $totgaps = 0;
	for (my $i=0; $i < scalar @seq; $i++) {
		if ($seq[$i] eq '-') {
            		$totgaps += 1;
			if ($in) {
				next;
			} else {
				push @start, $i;
				$in = 1;
			}
		} else {
			if ($in) {
				push @stop, $i-1;
				$in = 0;
			}
		}
	}
	## Account for sequences that end in a gap character!!! At this point they are missing a final stop value.
	if (scalar @start == (scalar @stop + 1) ) {
		push @stop , (scalar @seq - 1);
	}

	my %eliminate;

	if (scalar @start) {
		# reject those that are near other gaps
		if (scalar @start > 1) {
			for (my $i=0; $i<(scalar @start) - 1 ; $i++) {
				if ( ($start[$i+1] - $stop[$i]) < $proximity_th ) {
					$eliminate{$i} = 1;
					$eliminate{$i+1} = 1;
				}			
			}
		}
	
		# reject those greater than $maxDel threshold
		# reject those that are multiples of 3 with modulo function
		# reject those that extend to the beginning or end of the sequence
                # reject all if total gap length is greater than $maxPropExon of exon length
                
		for (my $i=0; $i<scalar @start; $i++) {
			my $length = $stop[$i] - $start[$i] + 1;
			$eliminate{$i} = 1 if $length > $maxDel;
			$eliminate{$i} = 1 unless ($length % 3);
			$eliminate{$i} = 1 if $start[$i] == 0; # Gap at beginning
			$eliminate{$i} = 1 if $stop[$i] == (scalar @seq - 1); # Gap at the end
            		$eliminate{$i} = 1 if $totgaps / scalar @seq >= $maxPropExon;
		}
		
		# Eliminate rejected gaps
		foreach my $i (sort {$b <=> $a} keys %eliminate) {
			splice @start , $i , 1;
			splice @stop , $i , 1;
		}
	}
	
	# return indepdent frameshift numbers
	return (scalar @start , \@start , \@stop);
}

