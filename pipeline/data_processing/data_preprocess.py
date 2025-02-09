import pandas as pd
import time
import argparse
from tqdm import tqdm, trange

parser = argparse.ArgumentParser()
parser.add_argument("--inputname", help="type input filename", type=str)
parser.add_argument("--trans", help="type trans data", type=str)
args = parser.parse_args()
inputname = args.inputname
trans = args.trans
mrna_275 = pd.read_csv(trans) 
mrna_275 = mrna_275[['Gene name', 'sequence']]
data = pd.read_csv('tmp/{}'.format(inputname), low_memory=False)
data = pd.merge(data, mrna_275, left_on='transcript_name', right_on='Gene name', how='inner')

try:
    data = data[['hybrid_seq', 'transcript_name', 'regulator_name', 'rem_tran_target_pos', 'remain_pos', 'on_reg_pos',
             'reg_hyb_target_pos', 'remain_seq', 'regulator_seq','pirscan_target_endpos', 'targeting_score',
             'raw_regulator_seq', 'idx', 'read_count', 'hybrid0', 'D', 'M', "pirscan binding site", "pirscan Target RNA sequence",
             'count', 'nor_readcount', 'nor_count', 'overlap','mir_init_pos', 'mir_end_pos','mir_energy',
             'mir_score', 'mir_target_pos', 'mir_transcript_seq', 'mir_regulator_seq','up_init_pos', 'up_end_pos',
             'RNAup_regulator_seq', 'RNAup_transcript_seq', 'RNAup_target_pos', 'RNAup_score', 'sequence']]
except:
    print('no pirscan1')
    data = data[['hybrid_seq', 'transcript_name', 'regulator_name', 'rem_tran_target_pos', 'remain_pos', 'on_reg_pos',
           'reg_hyb_target_pos', 'remain_seq', 'regulator_seq', 'idx', 'read_count', 'hybrid0', 'D', 'M',
           'count', 'nor_readcount', 'nor_count', 'overlap','mir_init_pos', 'mir_end_pos','mir_energy',
           'mir_score', 'mir_target_pos', 'mir_transcript_seq', 'mir_regulator_seq','up_init_pos', 'up_end_pos',
           'RNAup_regulator_seq', 'RNAup_transcript_seq', 'RNAup_target_pos', 'RNAup_score', 'sequence']]

data.reset_index(drop=True, inplace=True)
tmp_list = []
tmp2_list = []
tmp3_list = []
tmp4_list = []
tmp5_list = []
tmp6_list = []
tmp7_list = []
print(len(data))
time.sleep(1)
D_list = []
M_list = [] # regulator先沒加


for i in trange(len(data)):

    seq = data['sequence'][i]
    # pirscan
    try:
        pos_stop = int(data['pirscan_target_endpos'][i])+1
        pos_start = pos_stop-(len(data['regulator_seq'][i])-1)
        tmp_list.append('{}-{}'.format(pos_start, pos_stop))
        tmp2_list.append(seq[pos_start-1: pos_stop])
    except:
        pass 
    tmp3_list.append('{}_{}'.format(data['transcript_name'][i], data['regulator_name'][i]))
    
    # miranda
    mir_pos = data['mir_target_pos'][i].split('-')
    tmp4_list.append(seq[data['mir_init_pos'][i]+1+int(mir_pos[0])-1-1: data['mir_init_pos'][i]+1+int(mir_pos[1])-1])
    tmp5_list.append('{}-{}'.format(data['mir_init_pos'][i]+1+int(mir_pos[0])-1, data['mir_init_pos'][i]+1+int(mir_pos[1])-1)) 
    
    # RNAup
    up_pos = data['RNAup_target_pos'][i].split('-')
    tmp6_list.append(seq[data['up_init_pos'][i]+1+int(up_pos[0])-1-1: data['up_init_pos'][i]+1+int(up_pos[1])-1]) 
    tmp7_list.append('{}-{}'.format(data['up_init_pos'][i]+1+int(up_pos[0])-1, data['up_init_pos'][i]+1+int(up_pos[1])-1))
    
    # mutation是否在 overlap中
    D_raw = data['D'][i]
    M_raw = data['M'][i]
    if data['overlap'][i] == '0':
        D_list.append(D_raw)
        M_list.append(M_raw)

    else:
        if D_raw == '[]':
            D_list.append(D_raw)
        else:
            overlap = data['overlap'][i].split('-')
            hybrid_m = data['remain_pos'][i].split('-')
            on_m = data['rem_tran_target_pos'][i].split('-')
            over_start_m = int(overlap[0]) - int(hybrid_m[0]) + int(on_m[0])
            over_stop_m = int(overlap[1]) - int(hybrid_m[0]) + int(on_m[0])
            D = []
                    
            for m in eval(D_raw):
                if m < over_start_m or m > over_stop_m:
                    D.append(str(m))
                else:continue
            D_list.append(str(D))  
        if M_raw == '[]':
            M_list.append(M_raw)
        else:
            overlap = data['overlap'][i].split('-')
            hybrid_m = data['remain_pos'][i].split('-')
            on_m = data['rem_tran_target_pos'][i].split('-')
            over_start_m = int(overlap[0]) - int(hybrid_m[0]) + int(on_m[0])
            over_stop_m = int(overlap[1]) - int(hybrid_m[0]) + int(on_m[0])
            M = []
            for m in eval(M_raw):
                if m < over_start_m or m > over_stop_m:
                    M.append(str(m))
                else:continue
            M_list.append(str(M))   
            
print(len(D_list), len(M_list))        
try:
    data['pirscan_target_pos'] = tmp_list
    data['pir_target_mRNA_region'] = tmp2_list
except:
    print('no pirscan3')
data['mir_target_pos'] = tmp5_list
#data['mir_target_mRNA_region'] = tmp4_list
data['RNAup_target_pos'] = tmp7_list
#data['up_target_mRNA_region'] = tmp6_list

data['mRNA_len'] = [len(n) for n in data['sequence']]
data['hybrid_read'] = tmp3_list
data['D'] = D_list
data['M'] = M_list
data['A'] = [eval(str(data['D'][i]))+eval(str(data['M'][i])) for i in range(len(data))]
#data['mRNA_len'] = [len(n) for n in data['sequence']]
data = data.drop(columns=['sequence'])
outputname = inputname.split('/')[-1].replace('_detail_with_overlap.csv', '')
print(outputname)
data.to_csv('after_preprocess/{}_final.csv'.format(outputname), index=False)
