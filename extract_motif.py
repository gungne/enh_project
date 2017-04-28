
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
	pwm = pwm.astype(float)
	pfm = pfm.astype(float)
	# print(pfm)
	# print(np.shape(pfm[0,:])[1])
	for pos in range(0,np.shape(pfm[0,:])[1]):
		for nt in range(0,4):
			# print(pfm[nt,pos]/sum(pfm[:,pos]))
			if pfm[nt,pos] == 0:
				pfm[nt,pos] = sum(pfm[:,pos])*0.01

	for pos in range(0,np.shape(pfm[0,:])[1]):
		for nt in range(0,4):
			n_total = sum(pfm[:,pos])
			pwm[nt,pos] = pfm[nt,pos]/n_total
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
		# print(max_infomation,max_position)
		if  temp_information>max_infomation:
			
			max_infomation = temp_information
			max_position = n
	pwm_short =	pwm[:,max_position:max_position+l_kmer]	
	return pwm_short

def tf_proc(database_content,entry,l_kmer):
	ref_dict = dict()
	ref_dict_short = dict()
	if type(entry) is list:
		pass
	else:
		entry=[entry]
	for entry_name in entry:
		print(entry_name)
		pfm_comp = re.findall('\>.*'+entry_name+'\nA.*\nC.*\nG.*\nT.*',database_content,re.M|re.I)[0]
		tf_name,tf_pfm = pfm_parser(entry_name,pfm_comp)
		ref_dict[tf_name] = pfm2pwm(tf_pfm) 
		short_tf_pwm = shannon_trimmer(pfm2pwm(tf_pfm),l_kmer) 
		ref_dict_short[tf_name] = short_tf_pwm
	return ref_dict_short

def tf_proc_parse(database_content,entry,l_kmer):
	ref_dict = dict()
	# ref_dict_short = dict()
	if type(entry) is list:
		pass
	else:
		entry=[entry]
	for entry_name in entry:
		print(entry_name)
		pfm_comp = re.findall('\>.*'+entry_name+'\nA.*\nC.*\nG.*\nT.*',database_content,re.M|re.I)[0]
		tf_name,tf_pfm = pfm_parser(entry_name,pfm_comp)
		ref_dict[tf_name] = tf_pfm 
		# short_tf_pwm = shannon_trimmer(pfm2pwm(tf_pfm),l_kmer) 
		# ref_dict_short[tf_name] = short_tf_pwm
	return ref_dict

def fasta_parser(file):
	fasta_dict={}
	records = re.findall(">.*\n[ATCGatcg\n]*",file.read(),re.M)
	for fasta in records:
		id_sq = re.sub("[>\n]","",re.findall(">.*\n",fasta,re.M)[0])
		fasta_dict[id_sq] = re.sub("[\n]","",re.sub(">.*\n","",fasta)) 

	return fasta_dict
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

def nt_rconvert(num):
	if num == 0:
		row = 'A'
	if num == 1:
		row = 'C'
	if num == 2:
		row = 'G'
	if num == 3:
		row = 'T'
	return row

def kmer_score(pwm,kmer):
	# print(np.shape(pwm))
	total_score=1
	if np.shape(pwm)[1] != len(kmer):
		print('error')
	else:
		for index,nt in enumerate(kmer):
			# print(pwm[nt_convert(nt),index])
			# print(nt_convert(nt),index)
			total_score = total_score * pwm[nt_convert(nt),index]*100
			# print(total_score)
	return total_score
def pfm_writer(pfm,kmer):
	for index,nt in enumerate(kmer):
		# print(nt_convert(nt),index)
		pfm[nt_convert(nt),index] = pfm[nt_convert(nt),index]+1
	return pfm

def output_pfm_dict(dict_pfm):
	for item in dict_pfm:
		print(item)
		print('>'+item)
		indent = []
		for col_index,col_content in enumerate(dict_pfm[item].T):
			for posi in col_content.tolist():
				indent.append(len(str(max(posi))))
		# print(indent)
		for row_index,row_content in enumerate(dict_pfm[item][:]):
			# print(row_letter)
			# print(row_content)
			# for freq in row_content:
				# print(freq)
			row_list= row_content.tolist()
			# print(row_list)
			row_letter = nt_rconvert(row_index)
			print(row_letter +'  [') ,
			for col_index,spaces in enumerate(indent):
				# print(spaces)
				# print(row_list[col_index])
				# print(spaces),
				# print(len(str(row_list[0][col_index]))),
				for i in range(0,spaces-len(str(row_list[0][col_index]))):
					# print(spaces-len(str(row_list[0][col_index])))
					if spaces-len(str(row_list[0][col_index])) ==0:
						continue
					else:	
						sys.stdout.write(' '),
				sys.stdout.write(str(row_list[0][col_index])),
				sys.stdout.write(' '),
			sys.stdout.write(']')
		print('')
	return ''


if __name__ == "__main__":
	entry_name = sys.argv[1]
	fasta_name = sys.argv[2] 

	JASPAR_database = "sample_database.txt"
	handle_database = open(JASPAR_database)
	database_content = handle_database.read()

	handle_fasta = open(fasta_name)


	ref_pfm = tf_proc_parse(database_content,entry_name,6)
	# print(database_content)
	# single_motif = 
	# print(ref_dict_short)
	fasta_dict = fasta_parser(handle_fasta)
	l_kmer = 6
	# print(ref_pwm_short)
	new_pfm = np.zeros([4,6])
	print(ref_pfm)
	print(output_pfm_dict(ref_pfm))

	for read_nmer in fasta_dict:
		kmers = [ fasta_dict[read_nmer][n:n+l_kmer] for n in range(0,len(fasta_dict[read_nmer])-l_kmer+1)]
		best_kmer = ''
		best_score = 0
		for kmer in kmers:
			temp_score = kmer_score(ref_pwm_short[entry_name],kmer)
			# print(temp_score)
			if temp_score > best_score:
			# whatif the score are tied?
				best_kmer = kmer
				best_score = temp_score
		new_pfm = pfm_writer(new_pfm,best_kmer)
		# print(best_score)
	# print(output_pfm_dict(ref_pwm_short[entry_name]))
	print(new_pfm)
		# if best_kmer =='':
		# 	print(fasta_dict[read_nmer])
		# else:
		# 	print(best_kmer)	

	

























