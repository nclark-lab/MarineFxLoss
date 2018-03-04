#!/usr/bin/perl -w
use strict;

my $usage = "run_enrichments.pl <rank_file> <#>";
die $usage unless @ARGV == 2;

my $rank_file = $ARGV[0];
my $N = $ARGV[1];
my $p_value = 0.05;
my $th_count = 3; #below this value categories are removed.

system "mkdir results" unless (-e "results");

#set up topN file
system "head -n$N $rank_file > topN";

#Run thru all annotation files
foreach my $ann (qw(GNFtissueAnnotations curated.mSigDB MGIphenoExtended canonicalPath.mSigDB goBP.mSigDB)) {
	system "./find_enrichment_genes_generic.py topN -a enrichment_db/$ann -p $p_value -b $rank_file -o .txt -r $th_count";
	system "mv topN.genesEN.txt results/$rank_file\_$N\_$p_value\_$ann.txt";
}
exit;
