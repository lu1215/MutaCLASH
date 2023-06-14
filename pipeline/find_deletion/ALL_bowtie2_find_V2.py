import pandas as pd
import os
import argparse
import sys
import time

parser = argparse.ArgumentParser()
parser.add_argument("--inputname", help="type input filename", type=str)
parser.add_argument("--trans", help="type trans data", type=str)
args = parser.parse_args()
trans_data = args.trans
inputname = args.inputname

start = time.time()

def run_bwt2_trans(idx, remain_seq, trans_name, rem_tran_target_pos): # 將 identified region建成 index再將 hybrid read上 mRNA片段對上，最後取出bowtie2結果
    tmp = {}
    idx_t = idx
    with open('tmp_fasta1/read_test{}.fasta'.format(str(idx)), 'w') as f:
        f.write('>{}\n'.format(idx))
        f.write('{}\n'.format(remain_seq))
    with open('tmp_fasta1/ref_test{}.fasta'.format(str(idx)), 'w') as f:
        f.write('>{}\n'.format(trans_name))
        f.write('{}\n'.format(mrna_dict[trans_name][int(rem_tran_target_pos.split('-')[0])-1:int(rem_tran_target_pos.split('-')[1])]))

    os.system('bowtie2-build --quiet tmp_fasta1/ref_test{}.fasta tmp1/ref_test{}_tmp'.format(str(idx), str(idx)))
    os.system('bowtie2 --quiet -L 8 --gbar 1 -k 2 -i S,1,0 --score-min L,-18,0 --np 0 --mp 2,2 --rdg 3,1 --rfg 3,1 -x tmp1/ref_test{}_tmp -f tmp_fasta1/read_test{}.fasta --no-unal --no-hd -S bwt2_output1/test{}.sam'.format(str(idx), str(idx), str(idx)))
    sam = pd.read_csv('bwt2_output1/test{}.sam'.format(str(idx)),sep='\t', usecols=[0,3,5,17,18], header=None, names=['idx', 'pos_m', 'CIGAR_m', 'md1', 'md2'])
    os.system('rm tmp1/ref_test{}_tmp* tmp_fasta1/read_test{}.fasta tmp_fasta1/ref_test{}.fasta bwt2_output1/test{}.sam'.format(str(idx), str(idx), str(idx), str(idx)))
    z = 0
    if len(sam) == 0:
        pos = 1
        cigar = '{}M'.format(str(len(remain_seq)))
        md = 'MD:Z:{}'.format(str(len(remain_seq)))
    else:
        for i, row2 in sam.iterrows(): # 結果會按得分降序排列
            if 'D' in row2['CIGAR_m'] and 'I' not in row2['CIGAR_m']:
                z += 1
                idx = row2['idx']
                pos = row2['pos_m']
                cigar = row2['CIGAR_m']
                if row2['md2'][:2] == 'MD':
                    md = row2['md2']
                elif row2['md1'][:2] == 'MD':
                    md = row2['md1']
            else:continue
            if z == 1:
                break
        if z == 0:
            idx = sam['idx'][0]
            pos = sam['pos_m'][0]
            cigar = sam['CIGAR_m'][0]
            if sam['md2'][0][:2] == 'MD':
                md = sam['md2'][0]
            elif sam['md1'][0][:2] == 'MD':
                md = sam['md1'][0]
    
    tmp[idx_t] = [idx_t, pos, cigar, md]
    return tmp

def run_bwt2_reg(idx, remain_seq, rem_reg_target_pos):
    tmp = {}
    idx_t = idx
    with open('tmp_fasta2/read_test{}.fasta'.format(str(idx)), 'w') as f:
        f.write('>{}\n'.format(idx))
        f.write('{}\n'.format(remain_seq))
    with open('tmp_fasta2/ref_test{}.fasta'.format(str(idx)), 'w') as f:
        f.write('>{}\n'.format(idx))
        f.write('{}\n'.format(rem_reg_target_pos))
    os.system('bowtie2-build -q tmp_fasta2/ref_test{}.fasta tmp2/ref_test{}_tmp'.format(str(idx), str(idx)))
    os.system('bowtie2 --quiet -L 8 --gbar 1 -k 2 -i S,1,0 --score-min L,-18,0 --np 0 --mp 2,2 --rdg 3,1 --rfg 3,1 -x tmp2/ref_test{}_tmp -f tmp_fasta2/read_test{}.fasta --no-unal --no-hd -S bwt2_output2/test{}.sam'.format(str(idx), str(idx), str(idx)))
    sam = pd.read_csv('bwt2_output2/test{}.sam'.format(str(idx)),sep='\t', usecols=[0,3,5,17,18], header=None, names=['idx', 'pos_m', 'CIGAR_m', 'md1', 'md2'])
    os.system('rm tmp2/ref_test{}_tmp* tmp_fasta2/read_test{}.fasta tmp_fasta2/ref_test{}.fasta bwt2_output2/test{}.sam'.format(str(idx), str(idx), str(idx), str(idx)))
    z = 0
    if len(sam) == 0:
        pos = 1
        cigar = '{}M'.format(str(len(remain_seq)))
        md = 'MD:Z:{}'.format(str(len(remain_seq)))

    else:
        for i , row2 in sam.iterrows():
            if 'D' in row2['CIGAR_m'] and 'I' not in row2['CIGAR_m']:
                z += 1
                idx = row2['CIGAR_m']
                pos = row2['pos_m']
                cigar = row2['CIGAR_m']
                if row2['md2'][:2] == 'MD':
                    md = row2['md2']
                elif row2['md1'][:2] == 'MD':
                    md = row2['md1']
            else:continue        
            if z == 1:
                break
        if z == 0:
            idx = sam['idx'][0]
            pos = sam['pos_m'][0]
            cigar = sam['CIGAR_m'][0]
            if sam['md2'][0][:2] == 'MD':
                md = sam['md2'][0]
            elif sam['md1'][0][:2] == 'MD':
                md = sam['md1'][0]
    tmp[idx_t] = [idx_t, pos, cigar, md]
    return tmp

data = pd.read_csv(inputname)
print(len(data))

data['idx'] = [i for i in range(len(data))]
mrna = pd.read_csv(trans_data)
mrna_dict = {}
for i in range(len(mrna)):
    mrna_dict.update({mrna['Gene name'][i]: mrna['sequence'][i]})
     
import multiprocessing
from tqdm import trange

cores = multiprocessing.cpu_count()
print(cores-2)
pool = multiprocessing.Pool(processes=cores-2)
dic_trans = {}
dic_reg = {}
d_trans = []
d_reg = []
for i in trange(len(data)):
    idx = data['idx'][i]
    remain_seq = data['remain_seq'][i]
    trans_name = data['transcript_name'][i]
    rem_tran_target_pos = data['rem_tran_target_pos'][i]
    dicti_trans = pool.apply_async(run_bwt2_trans,(idx, remain_seq, trans_name, rem_tran_target_pos))
    d_trans.append(dicti_trans)
    
    #reg_on_hyb = data['reg_hyb_target_pos'][i].split('-')
    #remain_seq = data['hybrid_seq'][i][int(reg_on_hyb[0])-1: int(reg_on_hyb[1])]
    #reg_seq = data['regulator_seq'][i]
    #reg_pos = data['on_reg_pos'][i].split('-')
    #rem_reg_target_pos = reg_seq[int(reg_pos[0])-1: int(reg_pos[1])]
    #dicti_reg = pool.apply_async(run_bwt2_reg,(idx, remain_seq, rem_reg_target_pos))
    #d_reg.append(dicti_reg)

pool.close()
pool.join()

print('===== collectinf result =====')
print(len(d_trans))
for k in range(len(d_trans)):
    item = d_trans[k].get()
    dic_trans.update(item)
    #item = d_reg[k].get()
    #dic_reg.update(item)

data_bwt2_trans = pd.DataFrame.from_dict(dic_trans, orient='index', columns=['idx', 'pos_m', 'CIGAR_m', 'MD_m'])
#data_bwt2_reg = pd.DataFrame.from_dict(dic_reg, orient='index', columns=['idx', 'pos_r', 'CIGAR_r', 'MD_r'])

print(len(data_bwt2_trans))
data_bwt2_trans = data_bwt2_trans.reset_index(drop=True)
#data_bwt2_reg = data_bwt2_reg.reset_index(drop=True)
data = pd.merge(data, data_bwt2_trans, on='idx', how='inner')
#data = pd.merge(data, data_bwt2_reg, on='idx', how='inner')
print(len(data))
data.to_csv('ALL_output/{}.csv'.format(inputname.split('/')[-1].replace('.csv', '_step1')), index=False)
stop = time.time()
print('{} m'.format((stop - start)/60))
