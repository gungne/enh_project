
#==============================================================================
# batch download fastq base on whether two motifs are in the library 
#==============================================================================

import os
import sys
import re
import numpy as np
import Bio 
from subprocess import call


cur_dir = os.getcwd()
JASPAR_database = "pfm_vertebrates.txt"
exp_summary = "ERP008935_info.csv"

handle_expsum = open(exp_summary)
handle_database= open(JASPAR_database)

expsum_content = handle_expsum.read()

exp_list = re.findall('(?<=\n)ERR9[0-9]*',expsum_content )

# print(exp_list)
database_content = handle_database.read()

for exp_run in exp_list:
	# print(exp_run)
	run_detail = re.findall(str(exp_run)+'.*?(?<=,).*?40N.*?,',expsum_content )
	if len(run_detail)==0:
		continue
	run_name = re.findall('(?<=\,)[A-Za-z-0-9_]*?40N.*?(?=\,)',str(run_detail))
	# print(run_name)
	first_comp = re.findall('[A-Za-z0-9]{1}.*?(?=\_)',str(run_name))[0]
	second_comp = re.findall('[A-Za-z0-9]{1}.*?(?=\_)',str(run_name))[1]

	pfm_reg = '\>.*''\nA.*\nC.*\nG.*\nT.*'
	
	pfm_comp1 = re.findall('\>.*'+first_comp+'\nA.*\nC.*\nG.*\nT.*',database_content,re.M|re.I)
	pfm_comp2 = re.findall('\>.*'+second_comp+'\nA.*\nC.*\nG.*\nT.*',database_content,re.M|re.I)

	if len(pfm_comp1)>0 and len(pfm_comp2)>0:
		print(exp_run)
		call(['bash', 'sra_dl_single.sh', exp_run])
