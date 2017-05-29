# -*- coding: utf-8 -*-
"""
Created on Wed Mar 30 16:45:18 2017

@author: Gungnir
"""
#==============================================================================
# batch process the collapsed fasta file and assign with weights
#==============================================================================
import os
import sys
import re
import numpy as np

def fasta_parser(file):
	fasta_dict={}
	records = re.findall(">.*\n[ATCGatcg\n]*",file.read(),re.M)
	for fasta in records:
		id_sq = re.sub("[>\n]","",re.findall(">.*\n",fasta,re.M)[0])
		fasta_dict[id_sq] = re.sub("[\n]","",re.sub(">.*\n","",fasta)) 

	return fasta_dict

cur_dir = os.getcwd()
file_fa = sys.argv[1]
file_ref = sys.argv[2]

handle_fa = open(file_fa)
handle_ref = open(file_ref)


fa_dict = fasta_parser(handle_fa)
ref_content= handle_ref.read()
# ref_dict = dict()
for item in fa_dict:
	temp_item= re.sub('\|',',',item)
	temp_str =re.findall('n-R5s[0-9]*,'+temp_item,ref_content)[0]
	output = re.findall('n-R5s[0-9]*',temp_str)[0]
	print('>' + output)
	print(fa_dict[item])

 
	# for nmer in nmer_enriched_dict:
	# 	print('___',nmer)
	# 	for peak in nmer_enriched_dict[nmer]:
	# 		print('>')
	# 		print(peak)
	# print(nmer_enriched_dict)