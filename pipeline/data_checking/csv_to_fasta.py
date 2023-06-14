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
    data = data[['sequence', 'read_count']]
    with open(outputname, 'w') as f:
        for i in range(len(data)):
            f.writelines('>'+str(i)+'_'+str(data['read_count'][i])+'\n')
            f.writelines(data['sequence'][i]+'\n')

elif typename=='reference':
    data = data[['regulator_name', 'raw_regulator_seq']]
    with open(outputname, 'w') as f:
        for i in range(len(data)):
            f.writelines('>'+str(data['regulator_name'][i])+'\n')
            f.writelines(data['raw_regulator_seq'][i]+'\n')

#print('Program end with success.\n')