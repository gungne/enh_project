
#==============================================================================
# batch download fastq base on the 
#==============================================================================

import os
import sys
import re
import numpy as np
import Bio 
import math
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

def shannon_entropy(array_v):
	entropy = 0
	# print(array_v)
	for nt in array_v:
		# print(float(nt))
		if nt==0:
			nt=1
		entropy = entropy - float(nt)*math.log(nt,2)
	return entropy

def shannon_trimmer(pwm,l_kmer):
	shannon_pos=[]
	for pos in range(0,np.shape(pwm)[1]):	
		shannon_pos.append(shannon_entropy(pwm[:,pos]))
	max_infomation = 0   
	max_position =0
	# print(shannon_pos)
	for n in range(0,len(shannon_pos)-l_kmer):

		temp_information = 2*l_kmer -  sum(shannon_pos[n:n+l_kmer])

		print(max_infomation,max_position)
		if  temp_information>max_infomation:
			
			max_infomation = temp_information
			max_position = n
	pwm_short =	pwm[:,max_position:max_position+l_kmer]	
	return pwm_short

entry_name = sys.argv[1]

JASPAR_database = "pfm_vertebrates.txt"
handle_database = open(JASPAR_database)
database_content = handle_database.read()
ref_dict = dict()
ref_dict_short = dict()
pfm_comp = re.findall('\>.*'+entry_name+'\nA.*\nC.*\nG.*\nT.*',database_content,re.M|re.I)[0]
tf_name,tf_pfm = pfm_parser(entry_name,pfm_comp)
ref_dict[tf_name] = pfm2pwm(tf_pfm) 
short_tf_pwm = shannon_trimmer(pfm2pwm(tf_pfm),6) 
ref_dict_short[tf_name] = short_tf_pwm
print(ref_dict_short['PAX6'])
# single_motif = 
# print(pfm_comp)