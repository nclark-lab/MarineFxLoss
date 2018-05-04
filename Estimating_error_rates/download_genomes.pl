#!/usr/bin/perl

use strict;
open (C, "Chains");

my @chains=<C>;

my ($c, $genome);
my ($p1, $p2);
foreach $c (@chains){
    $c=~/mm10To(.*).over/;
    $genome=$1;
    $genome=~/([A-Z]{1}[a-z]*?)([A-Z].*)/;
    print "$genome\n";
    
    $genome=lc($1).$2;
print "1=$1\t2=$2\n";    
print "$genome\n";
    system "wget http://hgdownload.soe.ucsc.edu/goldenPath/$genome/bigZips/$genome.fa.gz";
}
