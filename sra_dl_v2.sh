#!/bin/bash
#=================================
# This script check if two components in the run list are in the the 
#=================================
run_dl_list=$(python pfm_exist_check.py)

for i in ${run_dl_list}
do
fastq-dump -A $i -O $HOME/enh_project/CAP_SELEX 
python find_repeat_exp.py $i CAP_SELEX/


done