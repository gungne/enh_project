#!/bin/bash
#=================================
# This script check if two components in the run list are in the the 
#=================================


dir="$HOME/enh_project/CAP_SELEX/" 
filename="$dir$1.fa"
if [ -f "$filename" ]
then
	echo "$filename found."
else
	fastq-dump -A $1 -O $HOME/enh_project/CAP_SELEX/ 
	echo "create $dir$1.fastq"
	# wc -l "$dir$1.fastq"
	fastq_to_fasta -v -Q33 -i "$dir$1.fastq" -o "$dir$1.fa"
	rm "$dir$1.fastq"
	# echo '$1.fa' >>'fasta_summary.txt'
	wc -l "$dir$1.fa" >> 'fasta_summary.txt'
fi

