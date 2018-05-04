#!/usr/bin/perl -w
use strict;
use PerlSubs;

# 58 mammals
my %species_include;
foreach (qw(hg19 panTro4 gorGor3 ponAbe2 nomLeu3 rheMac3 macFas5 papHam1 chlSab1 calJac3 saiBol1 otoGar3 tupChi1 speTri2 jacJac1 micOch1 criGri1 mesAur1 mm10 rn5 hetGla2 cavPor3 chiLan1 octDeg1 oryCun2 ochPri3 susScr3 vicPac2 camFer1 turTru2 orcOrc1 panHod1 bosTau7 oviAri3 capHir1 equCab2 cerSim1 felCat5 canFam3 musFur1 ailMel1 odoRosDiv1 lepWed1 pteAle1 pteVam1 myoDav1 myoLuc2 eptFus1 eriEur2 sorAra2 conCri1 loxAfr3 eleEdw1 triMan1 chrAsi1 echTel2 oryAfe1 dasNov3)) {
	$species_include{$_} = 1;
}

my %seq;
my ($g , $sp , $ex , $t_ex);
open (FILE , 'knownCanonical.exonNuc.45genes.nosep.fa');

my $first_time = 1;
while (<FILE>) {
	chomp;
	next if /^\s*$/; # skip pure whitespace lines

	# Break fasta header line into information components	
	if (/^>(\S+)_(\S+)_(\d+)_(\d+)\s/) {
		if ($first_time) {
			($g , $sp , $ex , $t_ex) = ($1 , $2 , $3 , $4);
			$first_time = 0;
		} else {
			($g , $sp , $ex , $t_ex) = ($1 , $2 , $3 , $4);
		}
	} else {
		if ($species_include{$sp}) {
			$seq{$g}->{$sp}->[$ex] = $_;
		}
	}
}
close FILE;
print STDERR "Finished reading in genes.\n";


open (FILE , 'knownCanonical.exonNuc.45genes.nosep.fa');
foreach my $line (read_input_clean("prestops.tsv")) {
	my @e = split /\t/ , $line;
	my $acc = pop @e;
	my ($g , $sp , $ex , $t_ex) = split /_/ , $acc;
	my $outline = join "\t" , @e , $seq{$g}->{$sp}->[$ex];
	print STDOUT "$outline\n";
}

print STDERR "Scans complete.\n";
close FILE;

exit;
