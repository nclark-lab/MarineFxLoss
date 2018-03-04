#!/usr/bin/perl
use strict;
use warnings;

my $startRun = time();

my $executable = '/home2/nclark/rap119/marine.bt.sim/BayesTraitsV3';
my $script1 = 'bt_indep';
my $script2 = 'bt_dep';

my $file_trait = 'trait_marine.txt';
#my $file_genes = '../pseudogenesUse.tsv';
#my $file_tree = 'tree_placentals.nex';
my $file_genes = $ARGV[0];
my $file_tree = $ARGV[1];


# Read and store gene presence/absence data in hash %data, store species names and order in @species
open( DATA, '<', $file_genes ) or die "Can't open $file_genes: $!";
my $firstLine = <DATA>;
my @species = split /\t/, $firstLine;
s{^\s+|\s+$}{}g foreach @species;
shift @species;	#remove "gene" column header
my $n_species = scalar @species;

my %data;
my %n1;
my %pattern;
my %tally;
foreach my $line (<DATA>) {
	chomp $line;
	my @array = split /\t/, $line;
	s{^\s+|\s+$}{}g foreach @array;
	my $gene = shift @array;
	my $n=0;

	my $string = join '' , @array;
	$data{$string} = \@array;

	foreach (@array) { $n++ if $_ == 1; }
	$n1{$string} = $n;

	$pattern{$gene} = $string;

	if (defined $tally{$string}) {
		$tally{$string}++;
	} else {
		$tally{$string} = 1;
	}
}
close DATA || die "Couldn't close file properly.";

print STDERR "Number of patterns: " , scalar keys %tally , "\n";
#foreach my $p (sort keys %tally) {
#	print STDERR "$p	$tally{$p}\n";
#}



# Read and store trait data in hash %trait
open( TRAIT, '<', $file_trait ) or die "Can't open $file_trait: $!";
my %trait;
foreach my $line (<TRAIT>) {
	chomp $line;
	my @array = split /\t/, $line;
	s{^\s+|\s+$}{}g foreach @array;
	$trait{$array[0]} = $array[1];
}
close TRAIT || die "Couldn't close file properly.";


my %result;
my $i=0;
foreach my $k (sort keys %data) {
	$i++;
	open( INFILE , '>', $file_genes.".infile" );
	for (my $i=0; $i<$n_species; $i++) {
		my $s = $species[$i];
		print INFILE "$s	$trait{$s}	$data{$k}->[$i]\n";
	}
	close INFILE;

	my @indep = `$executable $file_tree $file_genes."infile" < $script1`;
	my @dep   = `$executable $file_tree $file_genes."infile" < $script2`;
	chomp $indep[-2];
	$indep[-2] =~ s/\t$//;
	chomp $dep[-2];
	$dep[-2] =~ s/\t$//;
	
	$result{$k} = "$indep[-2]	$dep[-2]";
	print STDERR "$i ";
}
print STDERR "\n";

foreach my $g (sort keys %pattern) {
	my $p = $pattern{$g};
	print "$g	$n1{$p}	$result{$p}\n";		
}

my $endRun = time();
my $runTime = $endRun - $startRun;
print STDERR "\nRUNTIME	$runTime\n";

exit;
