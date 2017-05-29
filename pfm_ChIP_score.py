import os
import sys
import re
import numpy as np
import math
import extract_motif
from pfm_exist_check import find_component
from subprocess import call
from multiprocessing import Pool
import extract_spacing_pfm 
import collections
import matplotlib.pyplot as plt
def listdir_nohidden(path):
    for f in os.listdir(path):
        if not f.startswith('.'):
            yield f

def collapse_list_to_dict(rawlist):
	out_dict = dict()
	rawlist.sort()
	# print(rawlist)
	for item in rawlist:
		if item =='':
			continue
		if item in out_dict:
			out_dict[item] = out_dict[item] + 1
		else:
			out_dict[item] = 1
	return out_dict

def collapse_list(rawlist):
	rawlist.sort()
	newlist = []
	temp_item = ''
	for item in rawlist:
		if temp_item=='':
			temp_item = item
			newlist = newlist + [item]
			continue
		if item == temp_item:
			continue
		else:
			temp_item = item 
			newlist = newlist + [item]
	return newlist


if __name__ == "__main__":
	main_tf = sys.argv[1]
	fasta_handle = open(sys.argv[3])
	motif_sets_dir = sys.argv[2]

	JASPAR_database = "pfm_vertebrates.txt"
	handle_database = open(JASPAR_database)
	database_content = handle_database.read()

	ref_pfm = extract_motif.tf_proc_parse_s(database_content,main_tf,6)
	chip_data = extract_motif.fasta_parser(fasta_handle)
	file_list = listdir_nohidden(motif_sets_dir)

	pfm_dict = dict()

	for file_entry in file_list:
		file_entry_handle = open(motif_sets_dir + file_entry)
		file_entry_content = file_entry_handle.read()
		pfm_dict.update(extract_motif.tf_proc_all(file_entry_content))
	conc_chip_data = dict() 	

	broadpeaks = [] 
	# print(ref_pfm[main_tf],main_tf)
	for chip_entry in chip_data:
		broadpeaks = broadpeaks + extract_spacing_pfm.find_best_kmer_chip(ref_pfm[main_tf],chip_data[chip_entry])

	# collapse peaks
	# print(broadpeaks)
	chippeaks_dict = collapse_list_to_dict(broadpeaks)
	chippeaks_dict_ordered = collections.OrderedDict(sorted(chippeaks_dict.items()))
	# print(chippeaks_dict_ordered)







	enriched_list = []
	for collapsed_entry in chippeaks_dict_ordered:
		if chippeaks_dict_ordered[collapsed_entry] > 10:
			enriched_list  = enriched_list + [collapsed_entry]
	############plot distribution
	# histo_list = [] 
	# for collapsed_entry in chippeaks_dict_ordered:
	# 	histo_list = histo_list + [chippeaks_dict_ordered[collapsed_entry]]
	# plt.hist(histo_list,bins=max(histo_list),range= (0,max(histo_list)),log=True) 
	# plt.title("Chip_extracted peaks")
	# plt.xlabel("Occurence")
	# plt.ylabel("Frequency")
	# plt.show()
	############plot distribution
	chip_data = collections.OrderedDict(sorted(chip_data.items()))
	# print(enriched_list)
	# genome_region_list = []
	region_dict = dict() 
	for nmers in  enriched_list:
		if len(nmers)==26:
			for chip_entry in chip_data:
				temp_str = chip_data[chip_entry]
				temp_list = re.findall(nmers,temp_str)
				temp_count = 0 
				for item in temp_list:
					temp_count = temp_count + 1
					if nmers in region_dict: 
						region_dict[nmers] = region_dict[nmers] + [chip_entry]
					else:
						region_dict[nmers] = [chip_entry]
					# print(chip_entry,nmers,temp_count)
					# genome_region_list = genome_region_list + [chip_entry] 
				temp_list = re.findall(extract_spacing_pfm.reverse_comp(nmers),temp_str)
				temp_count = 0 
				for item in temp_list:
					temp_count = temp_count + 1
					if nmers in region_dict: 
						region_dict[nmers] = region_dict[nmers] + [chip_entry]
					else:
						region_dict[nmers] = [chip_entry]
					# print(chip_entry,nmers,temp_count)
					# genome_region_list = genome_region_list + [chip_entry] 


	# region_dict= collapse_list_to_dict(genome_region_list)
	for entry in region_dict:
		region_dict[entry] = collapse_list(region_dict[entry])

	nmer_enriched_dict = dict()
	for nmer in region_dict:
		nmer_enriched_dict[nmer]=[]
		for region in region_dict[nmer]:
			nmer_enriched_dict[nmer] = nmer_enriched_dict[nmer] + extract_spacing_pfm.find_best_kmer_chip(ref_pfm[main_tf],chip_data[region])
	# print(region_dict)
	# print(ref_pfm[main_tf])
	count= 0
	for nmer in nmer_enriched_dict:
		temp_list = []
		for peak in nmer_enriched_dict[nmer]:
			if len(peak) == 26:
				count=count+1
				temp_list = temp_list + [peak]
		nmer_enriched_dict[nmer]=temp_list
###########
	# for nmer in nmer_enriched_dict:
	# 	print('___',nmer)
	# 	for peak in nmer_enriched_dict[nmer]:
	# 		print('>')
	# 		print(peak)
	# print(nmer_enriched_dict)

###########
	enriched_pfm_dict = dict()
	for nmer in nmer_enriched_dict:
		for peak in nmer_enriched_dict[nmer]:
			enriched_pfm_dict = extract_spacing_pfm.pfm_dict_writer(enriched_pfm_dict, nmer, peak)
	# print(count)

	# print(enriched_pfm_dict)
	extract_motif.output_pfm_dict(enriched_pfm_dict,'ChIP_enriched.txt','output/' )











