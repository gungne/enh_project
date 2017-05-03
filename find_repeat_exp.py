import os
import sys
import re
import numpy as np
import Bio 
import math
import extract_motif
from pfm_exist_check import find_component
from subprocess import call
from extract_spacing_pfm import parse_pfm_dict
def reverse_comp(kmer):
	kmer_r=''
	for nt in kmer:
		if nt == 'A':
			kmer_r = kmer_r + 'T'
		if nt == 'C':
			kmer_r = kmer_r + 'G'
		if nt == 'G':
			kmer_r = kmer_r + 'C'
		if nt == 'T':
			kmer_r = kmer_r + 'A'
	return kmer_r[::-1]

def find_best_kmer(pwm,read_nmer,l_kmer):
	kmers = [ read_nmer[n:n+l_kmer] for n in range(0,len(read_nmer)-l_kmer+1)]
	best_kmer = ''
	best_score = 0
	for index,kmer in enumerate(kmers):
		temp_score = extract_motif.kmer_score(pwm,kmer)
		# print(temp_score)
		if temp_score > best_score:
		# whatif the score are tied?
			best_kmer_index = index
			best_kmer_orient = +1
			best_kmer = kmer
			best_score = temp_score
	for index,kmer in enumerate(kmers):
		temp_score = extract_motif.kmer_score(pwm,reverse_comp(kmer))
		# print(temp_score)
		if temp_score > best_score:
		# whatif the score are tied?
			best_kmer_index = index
			best_kmer_orient = -1
			best_kmer = reverse_comp(kmer)
			best_score = temp_score
	return best_kmer_index,best_kmer_orient,best_kmer,best_score

def nt_convert(char):
	if char == 'A':
		row = 0
	if char == 'C':
		row = 1
	if char == 'G':
		row = 2
	if char == 'T':
		row = 3
	return row

def pfm_dict_writer(pfm_dict, pfm_category, motif):
	if pfm_category in pfm_dict:
		pfm_dict[pfm_category]  = pfm_writer(pfm_dict[pfm_category],motif)
	else:
		pfm_dict[pfm_category] = np.zeros([4,len(motif)])
		# print(motif)
		# print(np.shape(motif))
		# print(pfm_dict[pfm_category])
		pfm_dict[pfm_category]  = pfm_writer(pfm_dict[pfm_category],motif)
	return pfm_dict

def pfm_writer(pfm,kmer):
	for index,nt in enumerate(kmer):
		# print(nt_convert(nt),index)
		pfm[nt_convert(nt),index] = pfm[nt_convert(nt),index]+1
	return pfm


		
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
for exp_entry in exp_entries:
	temp_exp = re.findall('ERR[0-9]*',exp_entry)[0]
	#check is the original selex data was there
	#check if there is some redundancy, if not generate ref_dict
	if os.path.isfile('Motif/'+temp_exp + '.fa'):
		continue
	else:
		if os.path.isfile(fasta_dir+temp_exp + '.fa'):
			continue
		else:
			call(['bash', 'sra_dl_single.sh', temp_exp])
		print(temp_exp)
		handle_fasta = open(fasta_dir+temp_exp + '.fa')
		fasta_dict = extract_motif.fasta_parser(handle_fasta)
		temp_dict = parse_pfm_dict(temp_exp,fasta_dir)
		extract_motif.output_pfm_dict(temp_dict,temp_exp+'.fa','Motif/')
		handle_fasta.close()
		call(['rm', fasta_dir+temp_exp+'.fa'])



handle_expsum.close()
handle_database.close()
