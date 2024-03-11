import pandas as pd
import argparse
import time
import re

parser = argparse.ArgumentParser()
parser.add_argument("--MD", help="type MD sam file", type=str)
parser.add_argument("--input", help="type input file", type=str)
parser.add_argument("--ref_col", help="type column of reference-name", type=str)
args = parser.parse_args()
MD_data = args.MD
input_data = args.input
ref_col = args.ref_col

start = time.time()

sam_dict = {
    # Mandatory Fields
    0: 'read_id',    # Query template NAME
    1: 'flag',       # bitwise FLAG
    2: 'ref_id',     # References sequence NAME
    3: 'target_pos', # 1-based leftmost mapping POSition
    4: 'score',      # MAPping Quality
    5: 'CIGAR',      # CIGAR string
    6: 'RNEXT',      # Ref. name of the mate/next read
    7: 'PNEXT',      # Position of the mate/next read
    8: 'TLEN',       # observed Template LENgth
    9: 'read_seq',   # segment SEQuence
    10: 'QUAL',      # ASCII of Phred-scaled base QUALity+33

    # Optional Fiels
    # Below is for <BWA>
    14: 'MD1',
    15: 'MD2',
}
cols = [0, 2, 3, 5, 9, 14, 15]

# read MD file
df = pd.read_csv(MD_data, sep='\t', header=None, usecols=cols, comment='@')
df = df.rename(columns=sam_dict)
df = df.rename(columns={'ref_id':ref_col, 'CIGAR':'CIGAR_m'})
df = df[df[ref_col]!='*']
df[['hybrid0','read_count']] = df['read_id'].str.split('_', expand=True)
df['hybrid0'] = df['hybrid0'].astype(int)

# find CIGAR tag
df['CIGAR_lst'] = df['CIGAR_m'].apply(lambda x: re.findall(r'[0-9]+[A-Z]',x))
for n in ['S','D','M','I']:
    df[n] = df['CIGAR_lst'].apply(lambda x: sum([ int(y[:-1]) for y in x if y.endswith(n)]))

df['CIGAR_m'] = df['CIGAR_m'].apply(lambda x: re.sub(r'[0-9]+[SH]', '', x))
df['remain_len'] = df['M'] + df['I']
df['end_pos'] = df['target_pos'] + df['read_seq'].str.len() + df['D'] - df['I'] - df['S'] - 1
df['rem_tran_target_pos'] = df['target_pos'].astype(str) + '-' + df['end_pos'].astype(str)
df['pos_m'] = 1

# find MD tag
df.loc[ df['MD1'].str[:4]=='MD:Z', 'MD_m'] = df['MD1']
df.loc[ df['MD2'].str[:4]=='MD:Z', 'MD_m'] = df['MD2']

del df['S']
del df['D']
del df['M']
del df['I']
del df['MD1']
del df['MD2']
del df['target_pos']
del df['read_id']
del df['read_count']
del df['end_pos']
del df['read_seq']
del df['CIGAR_lst']
df = df.drop_duplicates().reset_index(drop=True)

# read input file
data = pd.read_csv(input_data).drop_duplicates()
data['remain_len'] = data['remain_seq'].str.len()
data = pd.merge(data, df, on=['hybrid0','transcript_name','rem_tran_target_pos','remain_len']).drop_duplicates()
data['idx'] = data.index
del data['remain_len']

data.to_csv('ALL_output/{}.csv'.format(input_data.split('/')[-1].replace('.csv', '_step1')), index=False)
stop = time.time()
print('{} m'.format((stop - start)/60))
