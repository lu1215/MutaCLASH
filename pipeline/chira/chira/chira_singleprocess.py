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
    df = pd.DataFrame(list(zip(id, seq)), columns=['tagid', 'hybrid_seq'])
    return df

parser = argparse.ArgumentParser()
parser.add_argument("-i", help="type input filename", type=str)
parser.add_argument("-t", help="type target filename", type=str)
parser.add_argument("-s", help="singletons", type=str)
parser.add_argument("-o", help="type input outputname", type=str)
args = parser.parse_args()
inputname = args.i
targetname = args.t
sin = args.s
outputname = args.o

# read single and target file
target = read_fasta(targetname)
df_single = pd.read_csv(sin, sep='\t')
df_single.drop(columns=['geneid','symbol','region','groupid', 'tpm', 'tx_pos_strand','tx_pos_strand','length','genomic_pos','locus','score'], inplace=True)
df_single = df_single[df_single['txid'].isin(target['tagid'])].reset_index(drop=True)

# read seq data
df_seq = read_fasta(inputname)
data = pd.merge(df_seq, df_single, on='tagid')
data['transcript_name'] = data['txid']
data['tx_pos_start'] = data['tx_pos_start'] + 1
data['rem_tran_target_pos'] = data['tx_pos_start'].astype(str) + '-' + data['tx_pos_end'].astype(str)
data['init_pos'] = data['read_info'].str.split(',', expand=True)[0]
data['end_pos'] = data['read_info'].str.split(',', expand=True)[1]
data['remain_pos'] = data['init_pos'] + '-' + data['end_pos']
data['init_pos'] = data['init_pos'].astype(int)
data['end_pos'] = data['end_pos'].astype(int)
data['init_pos'] = data['init_pos'] - 1
data[['hybrid0','read_count']] = data['tagid'].str.split('_', expand=True)

# add "remain_seq"
lst = []
print(len(data))
for i, row in data.iterrows():
    print(i, end='\r')
    lst.append(row['hybrid_seq'][row['init_pos']: row['end_pos']])

data['remain_seq'] = lst
del data['tagid']
del data['txid']
del data['tx_pos_start']
del data['tx_pos_end']
del data['read_info']
del data['init_pos']
del data['end_pos']
data = data[['hybrid_seq','transcript_name','rem_tran_target_pos','remain_pos','remain_seq', 'read_count', 'hybrid0']]
data.to_csv(outputname, index=False)
