import argparse
import pandas as pd
import sys

parser = argparse.ArgumentParser()
parser.add_argument("--input", help="input filename", type=str)
parser.add_argument("--output", help="output filename", type=str)
parser.add_argument("--type", help="input or reference", type=str)
args = parser.parse_args()
inputname = args.input
outputname = args.output
typename = args.type

print('\n=================fasta_to_csv.py=================')
print('Input file:\n{}'.format(inputname))
print('Output file:\n{}'.format(outputname))
print('=================================================')

name = []
seq = []
col = []

with open(inputname, 'r') as f:
    lines = f.read().splitlines()
    for line in lines:
        if line[0]=='>':
            name.append(line[1:])
        else:
            seq.append(line)

if typename=='regulator':
    col = ['regulator_name','raw_regulator_seq']
elif typename=='transcript':
    col = ['Gene name','sequence']
else:
    print('[Error] Unknown type: {}'.format(typename))
    sys.exit(1)

data = pd.DataFrame(zip(name,seq), columns=col)
data = data.sort_values(by=col[0])
data.to_csv(outputname, index=False)

#print('Program end with success.\n')