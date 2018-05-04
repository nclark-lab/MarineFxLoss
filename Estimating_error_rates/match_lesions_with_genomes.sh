#!/bin/bash

./pseudo_caller_frameshift_batch_report.pl > shifts.tsv
./pseudo_caller_stops_batch_report.pl > prestops.tsv
./pull_nt.pl > stops.tsv
sort -t "\t" -k 1,3 stops.tsv shifts.tsv > sorted_list.tsv

./run_blast.pl

exit
