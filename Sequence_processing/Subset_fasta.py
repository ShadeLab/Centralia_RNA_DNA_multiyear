#!/usr/bin/python

import sys
import os
import pandas as pd
from Bio import SeqIO

# Arguments
print('Number of arguments: ', len(sys.argv), 'arguments.')
print('Input fasta: ', sys.argv[1])
print('List of sequences to keep: ', sys.argv[2])
print('Output fasta: ', sys.argv[3])
print('\n---\n')

# Read in sequences to keep
with open(sys.argv[2], 'r') as keepfile:
	seq_list = [i.replace('\n', '') for i in keepfile.readlines()]

# For each entry in the taxonomy fasta file, change header to be sintax approporiate
with open(sys.argv[3], 'w') as outfasta:
	for record in SeqIO.parse(sys.argv[1], 'fasta'):
		if record.id in seq_list:
			outfasta.write('>' + record.id + '\n')
			outfasta.write(str(record.seq) + '\n')
