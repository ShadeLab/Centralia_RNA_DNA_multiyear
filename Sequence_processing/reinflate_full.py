#!/usr/bin/python

import sys
import os
import pandas as pd
import numpy as np
from Bio import SeqIO

# Arguments
print('Number of arguments: ', len(sys.argv), 'arguments.')
print('Input fasta: ', sys.argv[1])
print('Input UC: ', sys.argv[2])
print('Output fasta: ', sys.argv[3])
print('\n---\n')

mapcolumns = ['type', 'cluster_n', 'length', 'percent_sim', 'orientation', 'none1', 'none2', 'none3', 'query', 'centroid']
mapping_file = pd.read_table(sys.argv[2], header=None, names=mapcolumns).drop(['type', 'cluster_n', 'length', 'percent_sim', 'orientation', 'none1', 'none2', 'none3'], axis=1).drop_duplicates()
mapping_file['centroid']=mapping_file['centroid'].apply(lambda x:str(x).replace('*','center'))
mapping_file['centroid'] = np.where(mapping_file['centroid'] == 'center', mapping_file['query'], mapping_file['centroid'])

seq_name_dict = {k: g["query"].tolist() for k,g in mapping_file.groupby("centroid")}

with open(sys.argv[3], 'w') as out_file:
	for fasta in SeqIO.parse(sys.argv[1], 'fasta'):
		unique_name, sequence = fasta.id, str(fasta.seq)
		unique_name = unique_name.split(';', 1)[0]
		for seq_name in seq_name_dict[unique_name]:
			sampleID = seq_name.rsplit('_', 1)[0]
			out_file.write('>' + seq_name + ';sample=' + sampleID + ';\n')
			out_file.write(sequence + '\n')
