#!/usr/bin/python

import sys
import os
import pandas as pd
from Bio import SeqIO

# Arguments
print('Number of arguments: ', len(sys.argv), 'arguments.')
print('Manifest file: ', sys.argv[1])
print('Output fastq: ', sys.argv[2])
print('\n---\n')

# Read in metadata
manifest_df = pd.read_csv(sys.argv[1], sep='\t')

# For each sample in manifest, read in the merged fastq file, change the read names to the sample name and count, then write to a new complete fastq file
with open(sys.argv[2], 'w') as outfastq:
	for index, row in manifest_df.iterrows():
		sample_id = row['SampleID'].replace('.', '_')
		print('Running ' + sample_id)
		count = 1
		for record in SeqIO.parse(row['merged_fastq'], 'fastq'):
			record.id = sample_id + '_' + str(count) + ' orig_name=' + record.id
			SeqIO.write(record, outfastq, 'fastq')
			count += 1
