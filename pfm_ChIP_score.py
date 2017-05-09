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

if __name__ == "__main__":

	fasta_handle = open(sys.argv[1])
	motif_sets_dir = sys.argv[2]
	