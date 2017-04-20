#!/bin/bash
for filename in $(eval echo "$1*.fastq")
do
	name=${filename##*/}
    base=${name%.fastq}
    echo "$1${base}.fastq"
	fastq_to_fasta -v -n -i "$1${base}.fastq" -o "$1${base}.fa"
	rm "$1${base}.fastq"
	# fastx_trimmer -v -f 1 -l 40 -i "$1${base}.fa" -o "$1${base}.trimmed.fa"
	# fastx_collapser -v -i "$1${base}.fa" -o "$1"collapsed/"${base}.collapsed.fa"
	# rm "$1${base}.fa"
	# head -n 1000 "$1"collapsed/"${base}.collapsed.fa" > "$1"enriched/"${base}.enriched.fa"
	# rm "$1"collapsed/"${base}.collapsed.fa"
	# meme <collapsed_ordered_list.fasta> -mod zoops -dna -minw 4 -maxw 20 -nmotifs 10 -maxsize 1000000.
done