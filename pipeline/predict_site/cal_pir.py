#-*- coding:utf-8 -*-
import os
import numpy as np
import pandas as pd
import sys
import time
import re
import multiprocessing
import argparse
#column sequence = mRNA_seq
def Reverse_complement(cal): # piRNA 序列 reverse complement
    rev=[]
    for i in range(0,len(cal)):
        rev.append(cal['raw_regulator_seq'][i][::-1]) ##revised1
    cal['pirev'] = rev
    ATCG_dict = {'A':'T','T':'A','C':'G','G':'C', 'U': 'A'}
    compl = []
    for j in range(0,len(rev)):
        revComSeqList = [ATCG_dict[k] for k in rev[j]]
        revComSeq = ''.join(revComSeqList)
        compl.append(revComSeq)
    cal['pirev_compl'] = compl    
    cal2 = cal.copy()
        
    return cal2
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

    data2=pd.DataFrame()
    data2['target_pos'] = data['rem_tran_target_pos']
    data2['pirevcom'] = data['pirev_compl']
    data2['mrna21'] = mrna21
    data2['mrnainitpos'] = initpos
    data2['mrnaendpos'] = endpos
    #data2['len_mrna21'] = [len(n) for n in mrna21]
    #data2['len_pirna'] = [len(n) for n in data['pirev_compl']]
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
            mrna21.append(row['sequence'][0:end_pos+expand_range+1]) #若往前取不到21個，則從前取到END+20
            initpos.append(0)
            endpos.append(end_pos+expand_range)
        elif len(row['sequence']) - end_pos < expand_range:
            mrna21.append(row['sequence'][init_pos-expand_range:])
            initpos.append(init_pos-expand_range)
            endpos.append(len(row['sequence'])-1)

    data2=pd.DataFrame()
    data2['target_pos'] = data['rem_tran_target_pos']
    data2['pirevcom'] = data['pirev_compl']
    data2['mrna21'] = mrna21
    data2['mrnainitpos'] = initpos
    data2['mrnaendpos'] = endpos
    #data2['len_mrna21'] = [len(n) for n in mrna21]
    #data2['len_pirna'] = [len(n) for n in data['pirev_compl']]
    return data2

def bestmatch(ori_mrna_initpos,ori_mrna_endpos,mrna,pirna,mrna_endpos,s): # pirscan
    # ori_mrna_endpos=copy.deepcopy(mrna_endpos)
    len_small_RNA = len(pirna)
    if len(mrna)>len_small_RNA-1:
        ret = []
        retpos = []
        for i in range(0,len(mrna)-(len_small_RNA-1)):
            ret.append(mrna[i:i+len_small_RNA])
            retpos.append(mrna_endpos)
            mrna_endpos = mrna_endpos -1 
            
        mismatch=[]
        onemis=[]
        score_list=[]
        score = 10
        gumis_zero_list = []
        nogumis_zero_list = []
        gumis_seedlist=[]
        nogumis_seedlist=[]
        gumis_list=[]
        nogumis_list=[]
        difer=[]
            
        gumis_zero = []
        nogumis_zero = []
        gumis_inseed=[]
        nogumis_inseed=[]
        gumis=[]
        nogumis=[]
        b={}
        
        
        pairing = []

        for j in range(0,len(ret)):
            for k in range(0,len_small_RNA):
                if pirna[k] != ret[j][k]: #mismatch
                    onemis.append(k+1)      #mismatch時先記錄位置
                    if k == 0:
                        if (pirna[k]=='C'and ret[j][k] =='T') or (pirna[k]=='A'and ret[j][k] =='G'):
                            gumis_zero.append(k+1)
                        
                        else:
                            nogumis_zero.append(k+1)
                            
                    elif k== 1 or k==2 or k==3 or k==4 or k==5 or k==6:
                        if (pirna[k]=='C'and ret[j][k] =='T') or (pirna[k]=='A'and ret[j][k] =='G'): #注意如果ret = G 對到的PI是A
                            gumis_inseed.append(k+1)                                                 #但ret=T 對到應該是C(complement是G)
                            score = score - 1.5
                        else:
                            nogumis_inseed.append(k+1)
                            score = score - 7
                    else:
                        if (pirna[k]=='C'and ret[j][k] =='T') or (pirna[k]=='A'and ret[j][k] =='G'): 
                            gumis.append(k+1)
                            score= score - 1.5
                        else:
                            nogumis.append(k+1)
                            score = score - 2
            mismatch.append(onemis)
            score_list.append(score)
            pairing.append(ret[j][:])
            gumis_zero_list.append(gumis_zero)
            nogumis_zero_list.append(nogumis_zero)
            gumis_seedlist.append(gumis_inseed)
            nogumis_seedlist.append(nogumis_inseed)
            gumis_list.append(gumis)
            nogumis_list.append(nogumis)
            
            onemis=[]
            score= 10
            gumis_zero=[]
            nogumis_zero=[]
            gumis_inseed=[]
            nogumis_inseed=[]
            gumis=[]
            nogumis=[]
        maxs =[ k for k in score_list if k == max(score_list)]    
        bestpos_s = [i for i,v in enumerate(score_list) if v == maxs[0]]
        if len(bestpos_s) != 1 : #multiple max score
            
            for m in range(0,len(bestpos_s)):
                tmpbestpos = bestpos_s[m]
                bestmrna_endpos = retpos[tmpbestpos]
                bestmrna_initpos = bestmrna_endpos - len_small_RNA-1 #len of mrna  # should flexible
                difer.append(min(ori_mrna_endpos,bestmrna_endpos)-max(ori_mrna_initpos,bestmrna_initpos)) #計算重疊部分
            bestpos_ss = [i for i,v in enumerate(difer) if v == max(difer)]
            
            bestpos = bestpos_s[bestpos_ss[-1]]    #if 與Target rna重疊部分一樣多，則保留位置前面的,找最前面的位置， 因為前面算分是從後開始算， 所以要取最後位
            bestmismatch = mismatch[bestpos]
            totalmismatch = len(bestmismatch)
            target_seq = pairing[bestpos]
                
            bestmrna_endpos = retpos[bestpos]
            # print(score_list)
            bestscore = max(score_list)
            # print(bestscore)
            gu_zero = gumis_zero_list[bestpos]
            nogu_zero = nogumis_zero_list[bestpos]

            gumis_seed=gumis_seedlist[bestpos]
            totalguseed =len(gumis_seed)
            nogumis_seed = nogumis_seedlist[bestpos]
            totalnguseed =len(nogumis_seed)

            gumis = gumis_list[bestpos]
            totalgu = len(gumis)
            nogumis = nogumis_list[bestpos]
            totalngu = len(nogumis)


            gu=np.append(gu_zero,gumis_seed)
            gu=np.append(gu,gumis)

            nogu=np.append(nogu_zero,nogumis_seed)
            nogu=np.append(nogu,nogumis)
                #s_i = str(s)+'_'+str(m)
            b[s] = [bestmrna_endpos,bestscore,target_seq,totalnguseed,totalguseed,totalngu,totalgu,totalmismatch,nogu,gu]

        else :    
            bestpos = bestpos_s[0]    
            bestmismatch = mismatch[bestpos]
            target_seq = pairing[bestpos]
            totalmismatch = len(bestmismatch)
        #  print(score_list)
            bestmrna_endpos = retpos[bestpos]
            bestscore = max(score_list)
        #   print(bestscore)
            gu_zero = gumis_zero_list[bestpos]
            nogu_zero = nogumis_zero_list[bestpos]

            gumis_seed=gumis_seedlist[bestpos]
            totalguseed =len(gumis_seed)
            nogumis_seed = nogumis_seedlist[bestpos]
            totalnguseed =len(nogumis_seed)

            gumis = gumis_list[bestpos]
            totalgu = len(gumis)
            nogumis = nogumis_list[bestpos]
            totalngu = len(nogumis)


            gu=np.append(gu_zero,gumis_seed)
            gu=np.append(gu,gumis)

            nogu=np.append(nogu_zero,nogumis_seed)
            nogu=np.append(nogu,nogumis)
            b[s] = [bestmrna_endpos,bestscore,target_seq,totalnguseed,totalguseed,totalngu,totalgu,totalmismatch,nogu,gu]
        return b
    else:
        #b = {}
        b = {}
        b[s] = [0,0,'',0,0,0,0,0,0,0]
        return b
if __name__ == "__main__":
    """

    input: step3_output.csv (have RNAup score cutoff)
           
    output: scoring function (only while regulator=piRNA)

    """
    start = time.time()
    parser = argparse.ArgumentParser()
    parser.add_argument("--inputname", help="type input filename", type=str)
    parser.add_argument("--ex", help="expand length", type=str)
    parser.add_argument("--reg", help="type reg data", type=str)
    parser.add_argument("--trans", help="type trans data", type=str)
   
    args = parser.parse_args()
    
    inputname = args.inputname
    reg_data = args.reg
    trans_data = args.trans
    ex_len = args.ex
    outputname = inputname.split('/')[-1].replace('.csv', '_scan.csv')
        
    parent_path = os.path.abspath(os.path.join(os.getcwd(), ".."))
    last_input_data = inputname#os.path.join(parent_path,inputname)
    output_path = '.'
    print('=====reading data=====')
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
    cal2 = Reverse_complement(lastdata_withseq)
    if ex_len == 'n':
        data2 = FindRet(cal2)
    else:
        data2 = FindRet_ex(cal2, int(ex_len))
    #data2.to_csv('data2.csv', index=False)    
    print('===== scoring =====')
    cores = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(processes=cores-2)# 維持執行的程序總數為processes，當一個程序執行完畢後會新增新的程序進去

    mrna_seq = data2['mrna21']
    pirna_seq = data2['pirevcom']
    mrna_pos = data2['mrnaendpos']
    b={}
    d=[]
    s= time.time()
    from tqdm import trange
    for j in trange(0,len(data2)):
        mrna = mrna_seq[j][::-1]  #因為比對計分PPT是從右到左，故為了方便 先反轉
        pirna = pirna_seq[j][::-1]
        mrna_endpos =mrna_pos[j]
        ori_mrna_endpos = int(cal2['rem_tran_target_pos'][j].split('-')[1])
        ori_mrna_initpos = int(cal2['rem_tran_target_pos'][j].split('-')[0])
        dicti = pool.apply_async(bestmatch,(ori_mrna_initpos,ori_mrna_endpos,mrna,pirna,mrna_endpos,j,))
        d.append(dicti)
    pool.close() # 關閉pool, 則不會有新的程序新增進去
    pool.join() # 必須在join之前close, 然後join等待pool中所有的執行緒執行完畢
    print('===== collecting result =====')
    for k in range(0,len(d)):
        items=d[k].get()
        b.update(items)
    data_score = pd.DataFrame.from_dict(b, orient='index')
    score_list = list(data_score[1])
    target_end_list = [str(n) for n in list(data_score[0])]
    target_seq = [str(n) for n in list(data_score[2])]
    #data_score = data_score.rename(columns={0:'target_score_endpos', 1:"targeting_score",2:'xgu_inseed',3:'gu_inseed',4:'xgu_innon-seed',5:'gu_innon-seed',6:'totalmismatch',7:'xGU_mispos',8: 'GU_mispos'})
    #DATA = pd.concat([lastdata, data_score], axis=1)
    #DATA= DATA.drop(columns=['pirev','pirev_compl','sequence'])
    lastdata['pirscan_target_endpos'] = target_end_list
    lastdata['targeting_score'] = score_list
    lastdata['-'] = '-'
    lastdata['pirscan Target RNA sequence'] = target_seq
    lastdata['pirscan_target_endpos'] = lastdata['pirscan_target_endpos'].astype('int')
    # lastdata['pirscan_target_initpos'] = lastdata['pirscan_target_endpos'] - len(lastdata['pirscan Target RNA sequence']) + 1
    lastdata['pirscan_target_initpos'] = lastdata['pirscan_target_endpos'] - lastdata['pirscan Target RNA sequence'].str.len() + 1
    lastdata['pirscan binding site'] = lastdata['pirscan_target_initpos'].astype('str') + lastdata['-'] + lastdata['pirscan_target_endpos'].astype('str')
    lastdata['pirscan score'] = score_list
    print('===== filtering score =====')
    #if score_cutoff is not None:
    #    DATA = DATA[DATA['targeting_score'].astype(float) <= score_cutoff]
    lastdata.to_csv('scan_output/'+outputname,index=False)
    print(time.time()-start)
    print(len(lastdata))

