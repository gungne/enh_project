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
file_tree = sys.argv[1]
file_ref = sys.argv[2]

handle_tree = open(file_tree)
handle_ref = open(file_ref)
weight_list = []

ori_tree= handle_tree.read()
new_tree = ori_tree
ref_fa = handle_ref.read()
reg_pattern_gi = "[A-Z]\w+[\.]?[0-9]?"
# gi_entry = re.search(reg_pattern_gi,line).group(0)
for gi_entry in re.findall(reg_pattern_gi, ori_tree):
	# print(gi_entry)
	reg_pattern_tag = ".*?\]"
	tag = re.search(gi_entry+reg_pattern_tag,ref_fa).group(0)
	short_tag = re.search("(?<=\[).*?(?=\])",tag).group(0)
	new_tree = re.sub(gi_entry,short_tag + "_" +gi_entry ,new_tree)

print(new_tree)