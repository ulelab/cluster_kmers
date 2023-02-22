#!/bin/bash
kmers=$(<test_kmers.txt)
peka_files=$(ls *5mer_distribution_UTR3.tsv)

python3 ../src/cluster_kmers.py \
-k $kmers \
-o ./results \
-p $peka_files \
-co seq_enrichment \
-n 3



