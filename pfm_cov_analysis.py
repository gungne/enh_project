import os
import sys
import re
import numpy as np
# import Bio 
import math
import extract_motif
from pfm_exist_check import find_component
from subprocess import call
from extract_spacing_pfm import parse_pfm_dict


def pfm_parser_v2(entry):
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

def tf_dict_parse(handle_fasta):
	ref_dict = dict()
	# ref_dict_short = dict()
	print(handle_fasta.read())
	for entry in handle_fasta.readlines()[::6]:
		print(entry)
		pfm_comp = re.findall('\>.*'+'\nA.*\nC.*\nG.*\nT.*',database_content,re.M|re.I)[0]
		tf_name,tf_pfm = pfm_parser(entry_name,pfm_comp)
		ref_dict[tf_name] = tf_pfm 
		# short_tf_pwm = shannon_trimmer(pfm2pwm(tf_pfm),l_kmer) 
		# ref_dict_short[tf_name] = short_tf_pwm
	return ref_dict

fasta_exprun = sys.argv[1] 
fasta_dir = sys.argv[2]
l_kmer = 6
JASPAR_database = "pfm_vertebrates.txt"
handle_database = open(JASPAR_database)
database_content = handle_database.read()

exp_summary = "ERP008935_info.csv"
handle_expsum = open(exp_summary)
expsum_content = handle_expsum.read()



# download bease on the existing fasta
first_comp,second_comp,pfm_comp1,pfm_comp2 = find_component(fasta_exprun,database_content,expsum_content)
# ERR[0-9]*?.*(?=ALX4_ALX4)
exp_entries=re.findall('ERR[0-9]*?.*'+first_comp+'_'+second_comp,expsum_content)


exp_dicts={}
for exp_entry in exp_entries:
	temp_exp = re.findall('ERR[0-9]*',exp_entry)[0]
	handle_fasta  = open('Motif/' + temp_exp +'.fa')
	# print(handle_fasta)
	pfm_dict = tf_dict_parse(handle_fasta)
	exp_dicts[temp_exp] = pfm_dict


print(exp_dicts)