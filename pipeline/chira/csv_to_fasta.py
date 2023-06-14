import argparse
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("--inputname", help="type input filename", type=str)
parser.add_argument("--outputname", help="type input filename", type=str)
args = parser.parse_args()
inputname = args.inputname
outputname = args.outputname

data = pd.read_csv(inputname)
data = data[['sequence', 'read_count']]

with open(outputname, 'w') as f:
    for i in range(len(data)):
        f.writelines('>'+str(i)+'_'+str(data['read_count'][i])+'\n')
        f.writelines(data['sequence'][i]+'\n')
