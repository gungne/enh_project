#!/bin/bash
for filename in $(eval echo "$1*.fastq")
do
	name=${filename##*/}
    base=${name%.fastq}
    echo "$1${base}.fastq"
    python find_pfm.py ${base} ERP008935_info.csv 

	# meme <collapsed_ordered_list.fasta> -mod zoops -dna -minw 4 -maxw 20 -nmotifs 10 -maxsize 1000000.
done