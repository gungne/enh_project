
#==============================================================================
# batch download fastq base on the 
#==============================================================================

import os
import sys
import re
import numpy as np
import Bio 
from subprocess import call

def pfm_parser(entry_name,entry):
	entry = str(entry)
	entry_fields = re.findall('(?<=\[)[0-9 ]*',entry)
	# print(entry_fields)
	for n in range(0,4):
		counts = [int(s) for s in entry_fields[n].split() if s.isdigit()]
		if n == 0:
			pfm=np.matrix(counts)
			# print(pfm)
		else:
			pfm=np.vstack([pfm,counts])
	return entry_name,pfm

def pfm2pwm(pfm):
	pwm = np.zeros(np.shape(pfm))
	pwm.astype(float)
	# print(np.shape(pfm[0,:])[1])
	for pos in range(0,np.shape(pfm[0,:])[1]):
		for nt in range(0,4):
			# print(pfm[nt,pos]/sum(pfm[:,pos]))
			pwm[nt,pos] = float(pfm[nt,pos])/float(sum(pfm[:,pos]))
	return pwm 

def shannon_trimmer(pwm,l_kmer):
	for pos in range(0,np.shape(pfm[0,:])[1]):






entry_name = sys.argv[1]

JASPAR_database = "pfm_vertebrates.txt"
handle_database = open(JASPAR_database)
database_content = handle_database.read()
ref_dict = dict()
pfm_comp = re.findall('\>.*'+entry_name+'\nA.*\nC.*\nG.*\nT.*',database_content,re.M|re.I)[0]
tf_name,tf_pfm = pfm_parser(entry_name,pfm_comp)
ref_dict[tf_name] = pfm2pwm(tf_pfm) 

print(ref_dict['PAX6'])
# single_motif = 
# print(pfm_comp)