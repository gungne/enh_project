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

cur_dir = os.getcwd()
file_run = sys.argv[1]
file_ref = sys.argv[2]
JASPAR_database = "pfm_vertebrates.txt"


handle_ref = open(file_ref)
handle_database= open(JASPAR_database)

run_detail = re.findall(file_run+'.*?(?<=,).*?40N.*?,',handle_ref.read())
run_name = re.findall('(?<=\,)[A-Za-z0-9_]*?40N.*?(?=\,)',str(run_detail))
first_comp = re.findall('[A-Za-z0-9]{1}.*?(?=\_)',str(run_name))[0]
second_comp = re.findall('[A-Za-z0-9]{1}.*?(?=\_)',str(run_name))[1]

# pfm_reg = '\>.*''\nA.*\nC.*\nG.*\nT.*'
# print('\>.*'+first_comp+'\nA.*\nC.*\nG.*\nT.*')
temp = handle_database.read()
#pfm_comp1 = re.findall('\>.*'+first_comp+'\nA.*\nC.*\nG.*\nT.*',handle_database.read(),re.M|re.I)

#pfm_comp2 = re.findall('\>.*'+first_comp+'\nA.*\nC.*\nG.*\nT.*',handle_database.read(),re.M|re.I)
pfm_comp1 = re.findall('\>.*'+first_comp+'\nA.*\nC.*\nG.*\nT.*',temp,re.M|re.I)
pfm_comp2 = re.findall('\>.*'+second_comp+'\nA.*\nC.*\nG.*\nT.*',temp,re.M|re.I)
# print(len(pfm_comp1))
if len(pfm_comp1)>0 and len(pfm_comp2)>0:
	print(pfm_comp1,pfm_comp2) 
else
	print(pfm_comp2)