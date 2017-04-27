import os
import sys
import re
import numpy as np
import Bio 
import math
import extract_motif
from pfm_exist_check import find_component
from subprocess import call
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
def output_pfm_dict(dict_pfm):
	for item in dict_pfm:
		print('>'+item)
		print('A  ['+ [ number+'\t' for nt in dict_pfm[item]])
fasta_exprun = sys.argv[1] 
fasta_dir = sys.argv[2]
l_kmer = 6
JASPAR_database = "pfm_vertebrates.txt"
handle_database = open(JASPAR_database)
database_content = handle_database.read()

exp_summary = "ERP008935_info.csv"
handle_expsum = open(exp_summary)
expsum_content = handle_expsum.read()

handle_fasta=open(fasta_dir+fasta_exprun + '.fa')
fasta_dict = extract_motif.fasta_parser(handle_fasta)

first_comp,second_comp,pfm_comp1,pfm_comp2 = find_component(fasta_exprun,database_content,expsum_content)

if len(pfm_comp1)>0 and len(pfm_comp2)>0:
	print(first_comp,second_comp)
else:
	exit()

new_pfm_dict = dict()
ref_pwm_short = extract_motif.tf_proc(database_content,[first_comp,second_comp],l_kmer)
new_pfm =np.zeros(np.shape(ref_pwm_short[first_comp]))


for read_nmer in fasta_dict:
	category_name= ''
	spacing=0
	# kmers = [ fasta_dict[read_nmer][n:n+l_kmer] for n in range(0,len(fasta_dict[read_nmer])-l_kmer+1)]
	best_kmer = ''
	best_score = 0
	fbest_kmer_index,fbest_kmer_orient,fbest_kmer,fbest_score = find_best_kmer(ref_pwm_short[first_comp],fasta_dict[read_nmer],l_kmer)
	rbest_kmer_index,rbest_kmer_orient,rbest_kmer,rbest_score = find_best_kmer(ref_pwm_short[second_comp],fasta_dict[read_nmer],l_kmer)
	if abs(fbest_kmer_orient*fbest_kmer_index-rbest_kmer_orient*rbest_kmer_index)>5+l_kmer or abs(fbest_kmer_orient*fbest_kmer_index-rbest_kmer_orient*rbest_kmer_index)<3:
		continue
	else:
		# print(fbest_kmer_orient,rbest_kmer_orient,fbest_kmer_index,rbest_kmer_index)
		if fbest_kmer_orient == 1 and rbest_kmer_orient == 1:
			if fbest_kmer_index < rbest_kmer_index:
				spacing = rbest_kmer_index-fbest_kmer_index-l_kmer
				category_name = first_comp+'::'+second_comp+'_'+'1+2+'+ '_' + str(spacing)
				motif = fasta_dict[read_nmer][fbest_kmer_index:rbest_kmer_index+l_kmer] 
			else:
				spacing = fbest_kmer_index-rbest_kmer_index-l_kmer
				category_name = first_comp+'::'+second_comp+'_''2+1+'+ '_' + str(spacing)
				motif = fasta_dict[read_nmer][rbest_kmer_index:fbest_kmer_index+l_kmer] 

		if fbest_kmer_orient == -1 and rbest_kmer_orient == -1:
			if fbest_kmer_index < rbest_kmer_index:
				spacing = rbest_kmer_index-fbest_kmer_index-l_kmer
				category_name = first_comp+'::'+second_comp+'_'+'2+1+'+ '_' + str(spacing)
				motif = reverse_comp(fasta_dict[read_nmer][fbest_kmer_index:rbest_kmer_index+l_kmer]) 
			else:
				spacing = fbest_kmer_index-rbest_kmer_index-l_kmer
				category_name = first_comp+'::'+second_comp+'_'+'1+2+'+ '_' + str(spacing)
				motif = reverse_comp(fasta_dict[read_nmer][rbest_kmer_index:fbest_kmer_index+l_kmer]) 

		if fbest_kmer_orient == 1 and rbest_kmer_orient == -1:
			if fbest_kmer_index < rbest_kmer_index:
				spacing = rbest_kmer_index-fbest_kmer_index-l_kmer
				category_name = first_comp+'::'+second_comp+'_'+'1+2-'+ '_' + str(spacing)
				motif = fasta_dict[read_nmer][fbest_kmer_index:rbest_kmer_index+l_kmer] 
			else:
				spacing = fbest_kmer_index-rbest_kmer_index-l_kmer
				category_name = first_comp+'::'+second_comp+'_'+'2-1+'+ '_' + str(spacing)
				motif = fasta_dict[read_nmer][rbest_kmer_index:fbest_kmer_index+l_kmer] 

		if fbest_kmer_orient == -1 and rbest_kmer_orient == 1:
			if fbest_kmer_index > rbest_kmer_index:
				spacing = fbest_kmer_index-rbest_kmer_index-l_kmer
				category_name = first_comp+'::'+second_comp+'_'+'1+2-' + '_' + str(spacing)
				motif = reverse_comp(fasta_dict[read_nmer][rbest_kmer_index:fbest_kmer_index+l_kmer]) 
			else:
				spacing = rbest_kmer_index-fbest_kmer_index-l_kmer
				category_name = first_comp+'::'+second_comp+'_'+'2-1+' + '_' + str(spacing)
				motif = reverse_comp(fasta_dict[read_nmer][fbest_kmer_index:rbest_kmer_index+l_kmer]) 
		# print(motif,category_name)
		# print(category_name,spacing)		
		# new_pfm = extract_motif.pfm_writer(new_pfm,best_kmer)
		# print(best_score)
		new_pfm_dict = pfm_dict_writer(new_pfm_dict, category_name, motif)
# print(ref_pwm_short[first_comp])
print(new_pfm_dict)

