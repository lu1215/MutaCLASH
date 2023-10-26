import argparse
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("--input", help="input filename", type=str)
parser.add_argument("--output", help="output filename", type=str)
parser.add_argument("--type", help="input or reference", type=str)
args = parser.parse_args()
inputname = args.input
outputname = args.output
typename = args.type

print('\n=================csv_to_fasta.py=================')
print('Input file:\n{}'.format(inputname))
print('Output file:\n{}'.format(outputname))
print('=================================================')

data = pd.read_csv(inputname)

if typename=='input':
    data['id'] = '>' + data.index.astype(str) + '_' + data['read_count'].astype(str)

elif typename=='regulator':
    data['id'] = '>' + data['regulator_name']
    data['sequence'] = data['raw_regulator_seq']

elif typename=='transcript':
    data['id'] = '>' + data['Gene name']

data = data[['id','sequence']]
data.to_csv(outputname, sep='\n', index=False, header=False)

#print('Program end with success.\n')