import os
import sys
import re
import numpy as np
import Bio 
from subprocess import call
import random 

def nt_num(num):
	if num == 0:
		nt = 'A'
	if num == 1:
		nt = 'T'
	if num == 2:
		nt = 'C'
	if num == 3:
		nt = 'G'
	return(nt)

consensus = sys.argv[1]
run= sys.argv[2]

for n in range(1,int(run)):
	print('> %s' % str(n));
	rand_seq = ''
	rand_seq_len = 40 - len(consensus)
	for m in range(1,rand_seq_len):
		rand_seq = rand_seq + nt_num(random.randint(0,3))
	insert_pos = random.randint(0,rand_seq_len)
	output_seq = rand_seq[0:insert_pos] + consensus + rand_seq[insert_pos:]
	print(output_seq)