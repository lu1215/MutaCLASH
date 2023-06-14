import argparse
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("--input", help="input filename", type=str)
parser.add_argument("--output", help="output filename", type=str)
args = parser.parse_args()
inputname = args.input
outputname = args.output

print('\n=================fasta_to_csv.py=================')
print('Input file:\n{}'.format(inputname))
print('Output file:\n{}'.format(outputname))
print('=================================================')

name = []
seq = []
with open(inputname, 'r') as f:
    lines = f.read().splitlines()
    for line in lines:
        if line[0]=='>':
            name.append(line[1:])
        else:
            seq.append(line)

data = pd.DataFrame(zip(name,seq), columns=['regulator_name','raw_regulator_seq'])
data = data.sort_values(by='regulator_name')
data.to_csv(outputname, index=False)

#print('Program end with success.\n')