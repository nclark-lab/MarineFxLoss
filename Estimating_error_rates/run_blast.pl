#!/usr/bin/perl -w
use strict;
use PerlSubs;

unless (-e "output") { system "mkdir output"; }
unless (-e "output_1line") { system "mkdir output_1line"; }
unless (-e "results_genes") { system "mkdir results_genes"; }
unless (-e "results_species") { system "mkdir results_species"; }


# Run thru each line
foreach my $line (read_input_clean("sorted_list.tsv")) {
#foreach my $line (read_input_clean("test_list.tsv")) {
	my @e = split /\t/ , $line;
	
	#capture gene and species from fields 0&2
	my ($g , $species , $pos) = ($e[0] , $e[2] , $e[5]);
	print STDERR "Processing $g $species $pos\n";

	#pull last element as query
	my $q = pop @e;
	#format string for query, rm '-'
	$q =~ s/-//g;
	#write query to tmp file
	my $name = "${g}_${species}_${pos}";
	output_string(">queries/tmp_$name",">$name\n$q\n");
	
	#execute BLAST with output format 0
	system "blastn -db ../blastdb/$species -query queries/tmp_$name -out output/tmp_$name.out -evalue 0.005 -outfmt '0' -num_descriptions 1 -num_alignments 1";
	
	#capture output
	my $blout = join '' , read_input("output/tmp_$name.out");
	my $snippet;
	if ( $blout =~ /Query=(.+)\n\nLambda/s ) {
		$snippet = "\nQuery=" . $1;
	}

	#execute BLAST outputformat 6
	system "blastn -db ../blastdb/$species -query queries/tmp_$name -out output_1line/tmp_$name.out -evalue 0.005 -outfmt '6 sseqid sseq' -max_target_seqs 1";
	my $snippet_1line = join '' , read_input("output_1line/tmp_$name.out");


	open OUTG , ">>results_genes/$g.out";
	open OUTS , ">>results_species/$species.out";

	#print separator, original line, and 2 BLAST snippets
	my $entry = "_" x 100;
	$entry .= "\n$line\n$snippet\n$snippet_1line\n";
	print OUTG $entry;
	print OUTS $entry;
	
	close OUTG;
	close OUTS;
	
}
exit;
