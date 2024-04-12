# This program is for merging read-count in FATSA file

import argparse
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("--input", help="type input filename", type=str)
parser.add_argument("--output", help="type output filename", type=str)
args = parser.parse_args()
inputname = args.input
outputname = args.output

print('\n=================merge_fasta.py==================')
print('Input file:\n{}'.format(inputname))
print('Output file:\n{}'.format(outputname))
print('=================================================')

# reassemble fasta file with read count
print('\n> Start reading fasta file and calculating read counts...')
print('Please wait a few seconds.')
seq = []
with open(inputname, 'r') as f:
    lines = f.read().splitlines()
    i = 0
    for line in lines:
        if line[0]!='>':
            seq.append(line)
            i += 1
data = pd.DataFrame(seq, columns =['sequence'])
data = data.sort_values(by='sequence')
data = data.groupby(['sequence']).size().reset_index(name='read_count')
data['index'] = data.index.astype(str)
data['name'] = '>' + data['index'] + '_' + data['read_count'].astype(str)
data = data[['name','sequence']]
print('End of calculaitng.')

# output fasta file with sequence & read count
print('\n> Start writting fasta file...')
data.to_csv(outputname, index=False, header=False, sep='\n')
print('Program end with success.\n')
