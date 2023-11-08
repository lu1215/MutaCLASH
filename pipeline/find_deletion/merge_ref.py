import pandas as pd
import argparse

def read_fasta(path):
    i = -1
    id = []
    seq = []
    with open(path, 'r') as f:
        for line in f:
            if line[0] == '>':
                id.append(line[1:-1])
                seq.append('')
                i += 1
            else:
                seq[i] += line[:-1]
    df = pd.DataFrame(list(zip(id, seq)), columns=['id', 'seq'])
    df['id'] = '>' + df['id']
    return df

parser = argparse.ArgumentParser()
parser.add_argument("--reg", help="type regulator file", type=str)
parser.add_argument("--tar", help="type target file", type=str)
args = parser.parse_args()

reg_data = read_fasta(args.reg)
tar_data = read_fasta(args.tar)
merge_data = pd.concat([reg_data,tar_data])
merge_data.to_csv('reference.fa', sep='\n', index=False, header=False)
