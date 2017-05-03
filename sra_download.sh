#!/bin/bash
for i in  $(eval echo {$1..$2})
do
	var1="ERR"
	combined="$var1$i"
	echo $combined
	fastq-dump -A $combined -O $HOME/enh_project/CAP_SELEX 

done

#ERR929466 ERR929721