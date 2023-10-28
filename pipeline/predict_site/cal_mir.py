import os
import numpy as np
import pandas as pd
import sys
import re
import argparse
import multiprocessing
from tqdm import trange
import time

def FindRet(cal2): # identified region前後最小延伸
    data = cal2.copy()
    mrna21=[]
    initpos=[]
    endpos=[]
    #expand_range = ex_len

    for i,row in data.iterrows(): #not sequence
        init_pos = int(row['rem_tran_target_pos'].split('-')[0])-1
        end_pos = int(row['rem_tran_target_pos'].split('-')[1])-1
        target_pos_len = end_pos - init_pos + 1
        if target_pos_len < len(row['raw_regulator_seq']):
            expand_range = len(row['raw_regulator_seq']) - target_pos_len
            if init_pos >= expand_range and len(row['sequence']) - end_pos >= expand_range:
                mrna21.append(row['sequence'][init_pos-expand_range:end_pos+expand_range+1])
                initpos.append(init_pos-expand_range)
                endpos.append(end_pos+expand_range)
            elif init_pos < expand_range and len(row['sequence']) - end_pos < expand_range:
                mrna21.append(row['sequence'][0:])
                initpos.append(0)
                endpos.append(len(row['sequence'])-1)
            elif init_pos < expand_range:
                mrna21.append(row['sequence'][0:end_pos+expand_range+1]) #若往前取不到21個，則從前取到END+20
                initpos.append(0)
                endpos.append(end_pos+expand_range)
            elif len(row['sequence']) - end_pos < expand_range:
                mrna21.append(row['sequence'][init_pos-expand_range:])
                initpos.append(init_pos-expand_range)
                endpos.append(len(row['sequence'])-1)
        else:
            mrna21.append(row['sequence'][init_pos:end_pos+1])
            initpos.append(init_pos)
            endpos.append(end_pos)
    
    
    data['mrna21'] = mrna21
    data['mir_init_pos'] = initpos
    data['mir_end_pos'] = endpos
    return data

def FindRet_ex(cal2, ex_len): # identified region前後固定延伸
    data = cal2.copy()
    mrna21=[]
    initpos=[]
    endpos=[]
    expand_range = ex_len

    for i,row in data.iterrows(): #not sequence
        init_pos = int(row['rem_tran_target_pos'].split('-')[0])-1
        end_pos = int(row['rem_tran_target_pos'].split('-')[1])-1
        if init_pos >= expand_range and len(row['sequence']) - end_pos >= expand_range:
            mrna21.append(row['sequence'][init_pos-expand_range:end_pos+expand_range+1])
            initpos.append(init_pos-expand_range)
            endpos.append(end_pos+expand_range)
        elif init_pos < expand_range and len(row['sequence']) - end_pos < expand_range:
            mrna21.append(row['sequence'][0:])
            initpos.append(0)
            endpos.append(len(row['sequence'])-1)
        elif init_pos < expand_range:
            mrna21.append(row['sequence'][0:end_pos+expand_range+1]) #若往前取不到21個，則從前取到END+20
            initpos.append(0)
            endpos.append(end_pos+expand_range)
        elif len(row['sequence']) - end_pos < expand_range:
            mrna21.append(row['sequence'][init_pos-expand_range:])
            initpos.append(init_pos-expand_range)
            endpos.append(len(row['sequence'])-1)
    
    
    data['mrna21'] = mrna21
    data['mir_init_pos'] = initpos
    data['mir_end_pos'] = endpos
    return data

def run_mir(idx, raw_regulator_seq, mrna21): # run miranda
    tmp = {}
    no = str(idx)
    with open('tmp_fasta/miRNA{}.fasta'.format(str(no)), 'w') as f:
        f.writelines('>'+no+'\n')
        f.writelines(raw_regulator_seq+'\n')
    with open('tmp_fasta/mRNA{}.fasta'.format(str(no)), 'w') as f:
        f.writelines('>'+no+'\n')
        f.writelines(mrna21+'\n')

    os.system('miranda tmp_fasta/miRNA{}.fasta tmp_fasta/mRNA{}.fasta -sc 0 -out mir_output/tmp_test/test{}.txt'.format(str(no), str(no), str(no)))
    os.system('grep \"Query\" mir_output/tmp_test/test{}.txt | grep \"3\'\" > mir_output/tmp_test/test{}_q.txt'.format(str(no), str(no)))
    os.system('grep \"Ref\" mir_output/tmp_test/test{}.txt | grep \"5\'\" > mir_output/tmp_test/test{}_r.txt'.format(str(no), str(no)))
    os.system('grep -A1 \"Scores for this hit:\" mir_output/tmp_test/test{}.txt | grep \">\" > mir_output/tmp_test/test_info{}.txt'.format(str(no), str(no)))

    data_info = pd.read_csv('mir_output/tmp_test/test_info{}.txt'.format(str(no)), delimiter = '\t', header=None, usecols=[0, 1, 2, 3, 5, 6], names=['regulator_name', 'transcript_name', 'Score', 'Energy', 'tran_pos', 'Align_len'])
    with open('mir_output/tmp_test/test{}_q.txt'.format(str(no)), 'r') as f:
        data_q = f.readlines()
    with open('mir_output/tmp_test/test{}_r.txt'.format(str(no)), 'r') as f:
        data_r = f.readlines()
    data_q = [item.replace(' ', '').replace('3', '').replace('5', '').split("'")[1] for item in data_q]
    data_r = [item.replace(' ', '').replace('3', '').replace('5', '').split("'")[1] for item in data_r]
    data_info['rem_tran_target_pos'] = data_info['tran_pos'].apply(lambda x:x.split(' ')[0]+'-'+x.split(' ')[1])
    data_info['regulator_name'] = data_info['regulator_name'].str.strip('>')
    data_info['transcript_seq'] = data_r
    data_info['regulator_seq'] = data_q
    os.system('rm tmp_fasta/miRNA{}.fasta tmp_fasta/mRNA{}.fasta mir_output/tmp_test/test{}.txt mir_output/tmp_test/test{}_q.txt mir_output/tmp_test/test{}_r.txt mir_output/tmp_test/test_info{}.txt'.format(no, no, no, no, no, no))
    if len(data_info) != 0: # miranda 有 data以最高分作為結果
        
        data_info = data_info.sort_values(by='Score')
        max_score = max(list(data_info['Score']))
        tmp_data = data_info[data_info['Score'] == max_score]
        if len(tmp_data) > 1: # 若最高分有多筆，以第一筆作為答案
            energy = (tmp_data['Energy'][0])
            score = (tmp_data['Score'][0])
            target = (tmp_data['rem_tran_target_pos'][0])
            transcript_seq = (tmp_data['transcript_seq'][0])
            regulator_seq = (tmp_data['regulator_seq'][0])

        else:
            energy = (tmp_data['Energy'][0])
            score = (tmp_data['Score'][0])
            target = (tmp_data['rem_tran_target_pos'][0])
            transcript_seq = (tmp_data['transcript_seq'][0])
            regulator_seq = (tmp_data['regulator_seq'][0])
    
    else:
        energy = ('-')
        score = ('-')
        target = ('-')
        transcript_seq = ('-')
        regulator_seq = ('-')
    tmp[idx] = [idx, energy, score, target, transcript_seq, regulator_seq]
    return tmp


parser = argparse.ArgumentParser()
parser.add_argument("--inputname", help="type input filename", type=str)
parser.add_argument("--reg", help="type reg data", type=str)
parser.add_argument("--trans", help="type trans data", type=str)
parser.add_argument("--ex", help="expand length", type=str)

args = parser.parse_args()
inputname = args.inputname
reg_data = args.reg
trans_data = args.trans
ex_len = args.ex

outputname = inputname.split('/')[-1].replace('.csv', '_mir.csv')
parent_path = os.path.abspath(os.path.join(os.getcwd(), ".."))
last_input_data = inputname#os.path.join(parent_path,inputname)
output_path = '.'
print(inputname)
print('=====reading data=====')

start = time.time()
lastdata = pd.read_csv(last_input_data)
mRNA = pd.read_csv(trans_data) #改mrna檔案
piRNA = pd.read_csv(reg_data) #改rna檔案
piRNA_dict = {}
for i in range(len(piRNA)):
    piRNA_dict.update({piRNA['regulator_name'][i]: piRNA['raw_regulator_seq'][i]})

print('===== mergeing sequence=====')
raw_seq_list = [piRNA_dict[n.split('_')[0]] for n in lastdata['regulator_name']]
lastdata['raw_regulator_seq'] = raw_seq_list
lastdata_withseq = pd.merge(lastdata,mRNA,left_on='transcript_name',right_on='Gene name',how='left')
print('===== finding complement and candidates =====')
if ex_len == 'n':   
    data2 = FindRet(lastdata_withseq) 
else:   
    data2 = FindRet_ex(lastdata_withseq, int(ex_len))
print('===== miRanda =====')

print(len(data2)) 
data2['idx'] = [i for i in range(len(data2))]
cores = multiprocessing.cpu_count()
pool = multiprocessing.Pool(processes=cores-2)

b = {}
d = []

for i in trange(len(data2)):
    idx = data2['idx'][i]
    raw_regulator_seq = data2['raw_regulator_seq'][i]
    mrna21 = data2['mrna21'][i]

    dicti = pool.apply_async(run_mir,(idx, raw_regulator_seq, mrna21))
    d.append(dicti)
pool.close()
pool.join()   

print('===== collecting result =====')
print(len(d))
for k in range(len(d)):
    items=d[k].get()
    b.update(items)

data_mir = pd.DataFrame.from_dict(b, orient='index', columns=['idx', 'mir_energy', 'mir_score', 'mir_target_pos', 'mir_transcript_seq', 'mir_regulator_seq'])     
data_mir = data_mir.reset_index()
data2 = data2.drop(columns=['sequence','mrna21'])
data2 = pd.merge(data2, data_mir, on='idx', how='inner')

data2.to_csv('mir_output/'+outputname, index=False)
stop = time.time()
print('{} m'.format((stop - start)/60))
