import pandas as pd
import argparse


parser = argparse.ArgumentParser()
parser.add_argument("-i1", help="type input filename1", type=str)
parser.add_argument("-i2", help="type input filename2", type=str)
parser.add_argument("-c", help="chimeras", type=str)
parser.add_argument("-o", help="type input outputname", type=str)

args = parser.parse_args()
inputname1 = args.i1
inputname2 = args.i2
chi = args.c
outputname = args.o


with open(inputname1, 'r') as f:
    wt1 = f.readlines()
          
    name_list = []
    seq_list = []

for i in range(0, len(wt1), 2):
    name_list.append(wt1[i][1:].replace('\n', ''))
    seq_list.append(wt1[i+1].replace('\n', ''))

wt1 = pd.DataFrame()
wt1['tagid'] = name_list
wt1['hybrid_seq'] = seq_list

del name_list, seq_list

with open(inputname2, 'r') as f:
    data = f.readlines()

data = [item.rstrip('\n') for item in data]

piRNA_dict = {}

for i in range(0, len(data), 2):
    piRNA_dict.setdefault(data[i].replace('>', ''), data[i+1])

data = pd.read_csv(chi, sep='\t')
data.drop(columns=['geneid1', 'geneid2', 'symbol1', 'symbol2','region1', 'region2', 'groupid1', 'groupid2',
                       'sequences', 'hybrid', 'hybrid_pos', 'mfe', 'tx_pos_strand1', 'tx_pos_strand2'], inplace=True)
data = pd.merge(wt1, data, on='tagid', how='inner')
print(len(data), len(list(set(data['hybrid_seq']))))

list1 = []
list2 = []
list3 = []
list4 = []
list5 = []
list6 = []
list7 = []
print(len(data))

for i, row in data.iterrows():
    print(i, end='\r')
    read_info = row['read_info'].split(',')
    list1.append('{}-{}'.format(row['tx_pos_start2']+1, row['tx_pos_end2']))
    list2.append('{}-{}'.format(read_info[2], read_info[3]))
    list3.append('{}-{}'.format(row['tx_pos_start1']+1, row['tx_pos_end1']))
    list4.append('{}-{}'.format(read_info[0], read_info[1]))
    list5.append(row['hybrid_seq'][int(read_info[2])-1: int(read_info[3])])
    list6.append(piRNA_dict[data['txid1'][i]])
                                        
data['rem_tran_target_pos'] = list1
data['remain_pos'] = list2
data['on_reg_pos'] = list3
data['reg_hyb_target_pos'] = list4

data['remain_seq'] = list5
data['regulator_seq'] = list6
data['read_count'] = [int(n.split('_')[1]) for n in data['tagid']] 
data['hybrid0'] = [int(n.split('_')[0]) for n in data['tagid']] 

data.drop(columns=['tx_pos_start1', 'tx_pos_end1', 'tx_pos_start2', 'tx_pos_end2'], inplace=True)
data = data.rename(columns={'txid1': 'regulator_name', 'txid2': 'transcript_name'})
data = data.drop(columns=['length1', 'length2','read_info', 'genomic_pos1', 'genomic_pos2', 'locus1', 'locus2', 'tpm1',
                                                       'tpm2', 'score1', 'score2', 'score','tagid'])
del list1, list2, list3, list4, list5, list6, list7
print(len(outputname))
data.to_csv(outputname, index=False)
