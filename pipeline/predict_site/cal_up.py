import pandas as pd
import os

import time
import argparse
import multiprocessing

def FindRet(cal2): # identified region前後最小延伸
    data = cal2.copy()
    mrna21=[]
    initpos=[]
    endpos=[]
    #expand_range = 20
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
                mrna21.append(row['sequence'][0:end_pos+expand_range+1]) 
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
    data2=pd.DataFrame()
    data2['target_pos'] = data['rem_tran_target_pos']
    data2['mrna21'] = mrna21
    data2['init_pos'] = initpos
    data2['end_pos'] = endpos
    return data2

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
            mrna21.append(row['sequence'][0:end_pos+expand_range+1]) 
            initpos.append(0)
            endpos.append(end_pos+expand_range)
        elif len(row['sequence']) - end_pos < expand_range:
            mrna21.append(row['sequence'][init_pos-expand_range:])
            initpos.append(init_pos-expand_range)
            endpos.append(len(row['sequence'])-1)

    data2=pd.DataFrame()
    data2['target_pos'] = data['rem_tran_target_pos']
    data2['mrna21'] = mrna21
    data2['init_pos'] = initpos
    data2['end_pos'] = endpos
    return data2

def RNAup(idx, trans_seq, reg_seq, target_region): # RNAup 指令與結果
    seq = trans_seq+'&'+reg_seq
    #print(seq)
    f = os.popen("echo '"+seq+"'|RNAup -b -d2 --noLP -c 'S' -o")
    up = f.readlines()
    out = up[0].split('  ')
    up_seq = up[1].replace('\n', '')
    #print(out)
    s = out[0]
    trans_pos = out[1].replace(' ', '').split(',')
    reg_pos = out[3].replace(' ', '').split(',')
    score = out[4].split('=')[0].replace(' ', '').replace('(', '')
    tmp = cal_pos(idx, s, trans_seq, reg_seq, trans_pos, reg_pos, target_region, score, up_seq)
    return tmp

def cal_pos(idx, s, trans_seq, reg_seq, trans_pos, reg_pos, target_region, score, up_seq): # binding site計算
    #print(idx, s, trans_seq, reg_seq, trans_pos, reg_pos, target_region, score)
    # cal pos
    tmp = {}
    no = str(idx)
    if up_seq == '&':
        tmp[idx] = [reg_seq, trans_seq, '0-0', 100]
    else:
        trans_s = s.split('&')[0][::-1].replace('(', ')')
        reg_s = s.split('&')[1]
        trans_bind_seq = trans_seq[int(trans_pos[0])-1:int(trans_pos[1])][::-1]
        reg_bind_seq = reg_seq[int(reg_pos[0])-1:int(reg_pos[1])]
        seq1_len = len(trans_bind_seq)
        i = 0
        while i <= seq1_len-1:
            #print(i, seq1_len)
            try:
                if trans_s[i] != reg_s[i]:
                    if trans_s[i] == '.':
                        reg_bind_seq = reg_bind_seq[:i]+'-'+reg_bind_seq[i:]
                        reg_s = reg_s[:i]+'.'+reg_s[i:]
                    else:
                        #print(trans_s[i], reg_s[i])
                        trans_bind_seq = trans_bind_seq[:i]+'-'+trans_bind_seq[i:]
                        trans_s = trans_s[:i]+'.'+trans_s[i:]
                        seq1_len += 1
                else:
                    i += 1
            except:
                 i += 1
        p_front = 0
        if int(reg_pos[0]) != 1:
            ex_len = int(reg_pos[0])-1
            for i in range(ex_len):
                reg_bind_seq = reg_seq[ex_len-i-1] + reg_bind_seq
                if len(trans_seq)-(i+1) >= int(trans_pos[1]):
                    trans_bind_seq = trans_seq[int(trans_pos[1])+(i)] + trans_bind_seq
                    p_front += 1
                else:
                    trans_bind_seq = '-' + trans_bind_seq

        p_back = 0
        if int(reg_pos[1]) != len(reg_seq):
            ex_len = len(reg_seq) - int(reg_pos[1])
            for i in range(ex_len):
                reg_bind_seq = reg_bind_seq + reg_seq[len(reg_seq)-(ex_len-i)]
                if int(trans_pos[0])-(i+1) >= 1:
                    trans_bind_seq = trans_bind_seq + trans_seq[int(trans_pos[0])-(i+1)-1]
                    p_back += 1
                else:
                    trans_bind_seq = trans_bind_seq + '-'
        target_pos = '{}-{}'.format(int(target_region.split('-')[0])+(int(trans_pos[0])-p_back)-1,
                                int(target_region.split('-')[1])-(len(trans_seq)-int(trans_pos[1])-p_front))

        tmp[idx] = [idx, reg_bind_seq[::-1], trans_bind_seq[::-1], target_pos, score]
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

outputname = inputname.split('/')[-1].replace('.csv', '_RNAup.csv')

print(inputname)
print('=====reading data=====')
start = time.time()
data = pd.read_csv(inputname)
mrna_275 = pd.read_csv(trans_data) 
mrna_275 = mrna_275[['Gene name', 'sequence']]
data = pd.merge(data, mrna_275, left_on='transcript_name', right_on='Gene name', how='inner')
print(len(data))
if ex_len == 'n':   
    data2 = FindRet(data) 
else:   
    data2 = FindRet_ex(data, int(ex_len))

list6 = []
for i in range(len(data)):
    list6.append('1-{}'.format(len(data2['mrna21'][i])))
    
data['up_init_pos'] = data2['init_pos']
data['up_end_pos'] = data2['end_pos']
data['rem_tran_target_pos2'] = list6
data['sequence'] = data2['mrna21']
data = data.rename(columns={'sequence': 'transcript_seq'})
del data2

print(len(data))
data['idx'] = [i for i in range(len(data))]
cores = multiprocessing.cpu_count()
pool = multiprocessing.Pool(processes=cores-2)
b = {}
d = []

for l in range(len(data)):
    #print(l)
    idx = data['idx'][l]
    trans_seq = data['transcript_seq'][l]
    reg_seq = data['regulator_seq'][l]
    target_region = data['rem_tran_target_pos2'][l]
    dicti = pool.apply_async(RNAup,(idx, trans_seq, reg_seq, target_region))
    #dicti = RNAup(idx, trans_seq, reg_seq, target_region)
    d.append(dicti)
pool.close()
pool.join()
#print(d)
print('===== collecting result =====')
print(len(d))
for k in range(len(d)):
    items=d[k].get()
    b.update(items)
data_up = pd.DataFrame.from_dict(b, orient='index', columns=['idx', 'RNAup_regulator_seq', 'RNAup_transcript_seq', 'RNAup_target_pos', 'RNAup_score'])
data_up = data_up.reset_index()
data = pd.merge(data, data_up, on='idx', how='inner')
data.to_csv('up_output/'+outputname, index=False)
stop = time.time()
print('{} m'.format((stop - start)/60))
