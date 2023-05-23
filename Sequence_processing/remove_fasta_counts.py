#!/usr/bin/python

import sys
import os
import pandas as pd
from Bio import SeqIO

# Arguments
print('Number of arguments: ', len(sys.argv), 'arguments.')
print('Input fasta: ', sys.argv[1])
print('Output fastq: ', sys.argv[2])
print('\n---\n')

# Remove the count value (size=) from the fasta headers
with open(sys.argv[2], 'w') as outfasta:
	for record in SeqIO.parse(sys.argv[1], 'fasta'):
		new_id = record.id.split(';')[0]
		outfasta.write('>' + new_id + '\n')
		outfasta.write(str(record.seq) + '\n')
