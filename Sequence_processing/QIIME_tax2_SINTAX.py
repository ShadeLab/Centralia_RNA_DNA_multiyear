#!/usr/bin/python

import sys
import os
import pandas as pd
from Bio import SeqIO

# Arguments
print('Number of arguments: ', len(sys.argv), 'arguments.')
print('Taxonomy fasta: ', sys.argv[1])
print('Taxonomy table: ', sys.argv[2])
print('Output taxonomy fasta: ', sys.argv[3])
print('\n---\n')

# Read in taxonomy table
tax_df = pd.read_csv(sys.argv[2], sep='\t', header=None)
tax_df.columns = ['Reference', 'Taxonomy']
tax_dict = {}
for index, row in tax_df.iterrows():
	taxonomy = row['Taxonomy'].replace('D_0__', 'tax=d:').replace(';D_1__', ',p:').replace(';D_2__', ',c:').replace(';D_3__', ',o:').replace(';D_4__', ',f:').replace(';D_5__', ',g:').replace(';D_6__', ',s:')
	tax_dict[row['Reference']] = taxonomy

# For each entry in the taxonomy fasta file, change header to be sintax approporiate
with open(sys.argv[3], 'w') as outfasta:
	for record in SeqIO.parse(sys.argv[1], 'fasta'):
		outfasta.write('>' + record.id + ';' + tax_dict[record.id] + ';\n')
		outfasta.write(str(record.seq) + '\n')
