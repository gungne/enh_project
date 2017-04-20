#!/bin/bash
for filename in $(eval echo "$1*.fa")
do
	name=${filename##*/}
    base=${name%.enriched.fa}
    # echo "$1${base}.weighted.fa"
    python weight_collapsed.py $1${base}.enriched.fa | cat - $1${base}.enriched.fa > $1weighted/${base}.weighted.fa

done