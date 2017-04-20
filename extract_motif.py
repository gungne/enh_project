
#==============================================================================
# batch download fastq base on the 
#==============================================================================

import os
import sys
import re
import numpy as np
import Bio 
from subprocess import call

def pfm_parser(entry):
	pfm = np.array([]) 
	entry = str(entry)
	entry_fields = re.findall('(?<=\[)[0-9 ]*',entry)
	
	for n in range(0,3):
		counts = re.findall("[0-9]*",entry_fields[n])
		for pos in counts:
			pfm[n:].append(pos)
	print(pfm)
# def shannon_trimmer():




entry_name = sys.argv[1]

JASPAR_database = "pfm_vertebrates.txt"
handle_database = open(JASPAR_database)
database_content = handle_database.read()

pfm_comp = re.findall('\>.*'+entry_name+'\nA.*\nC.*\nG.*\nT.*',database_content,re.M|re.I)[0]
pfm_parser(pfm_comp)
# single_motif = 
# print(pfm_comp)