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
file = sys.argv[1]
handle= open(file)

weight_list = []

for line in handle.readlines():
	if line.startswith('>'):
		reg_pattern = "(?<=-)[0-9]*"
		weight = re.search(reg_pattern,line).group(0)
		weight_list.append(int(weight))

print('>WEIGHTS', end=" ")
for item in weight_list: 
	print('%.4f ' % (item/weight_list[0]), end=" ")

print()
