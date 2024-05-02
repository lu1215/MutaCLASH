import pandas as pd
import sys
import numpy as np


def gu_fun(rem1,rem2,RNA_Seq,pirna_com,target_score,rem_seq,mRNA_start,mRNA_end):

    rem2=rem2+1
    multi=[]
    if rem2 < 0:
        return ["0-0",0,0,0,0,0,0,"0","0"]
    
                            
    for i in range(rem1,rem2):
        DY=RNA_Seq[i:i+21]
        DY=DY[::-1]
        DY=list(DY)
        gu=[]
        non_gu=[]
        match=[]


        for j in range(21):
            if [DY[j],pirna_com[j]]==['G','A'] or [DY[j],pirna_com[j]]==['T','C'] :
                gu.append(j+1)

            elif DY[j]==pirna_com[j] :
                match.append(j+1)
            else:
                non_gu.append(j+1)
        gu_inseed=sum(8>k>1 for k in gu)
        gu_outseed=sum(k>7 for k in gu)
        nongu_inseed=sum(8>k>1 for k in non_gu)
        nongu_outseed=sum(k>7 for k in non_gu)

        score= 10-7*nongu_inseed-1.5*gu_inseed-2*nongu_outseed-1.5*gu_outseed
        target_score.append(score)
    max_score=max(target_score,default=0)
    max_score_p=[i for i,v in enumerate(target_score) if v==max_score]

    local_mRNAstart=mRNA_start-rem1
    local_mRNAend=mRNA_end-rem1
    pos_dis=-1
    close_pos=0
    sta=0
    end=0
    for p in range(0,len(max_score_p)):
        sta=local_mRNAstart if max_score_p[p] < local_mRNAstart else max_score_p[p]
        end=local_mRNAend if max_score_p[p]+20 > local_mRNAend else max_score_p[p]+20

        if (end-sta) > pos_dis:
            pos_dis=end-sta
            close_pos=p

    max_score_p=max_score_p[close_pos:close_pos+1]

    for max_score_p in max_score_p:
        target_position='{}-{}'.format(max_score_p+rem1+1,max_score_p+rem1+21)
        rna_position='{}-{}'.format(rem_seq[0]+1,rem_seq[1])



        DYl=RNA_Seq[max_score_p+rem1:max_score_p+rem1+21]
        DYl=DYl[::-1]
        DYl=list(DYl)
        gul=[]
        non_gul=[]
        matchl=[]

        for i in range(21):
            if [DYl[i],pirna_com[i]]==['G','A'] or [DYl[i],pirna_com[i]]==['T','C'] :
                gul.append(i+1)

            elif DYl[i]==pirna_com[i] :
                matchl.append(i+1)
            else:
                non_gul.append(i+1)
        gu_inseedl=sum(8>i>1 for i in gul)
        gu_outseedl=sum(i>7 for i in gul)
        nongu_inseedl=sum(8>i>1 for i in non_gul)
        nongu_outseedl=sum(i>7 for i in non_gul)
        total_mismatch=len(gul)+len(non_gul)
        em=[]
        
        multi.append([target_position,max_score,nongu_inseedl,gu_inseedl,nongu_outseedl,gu_outseedl,total_mismatch,str(non_gul)[1:-1],str(gul)[1:-1]])

    return multi[0]



def complement(seq):
    RNA=0
    result_seq=[]
    if "U" in seq:
        RNA=1
    for i in seq:
        if RNA:
            # RNA
            if i=="G":
                result_seq.append("C")
            elif i=="C":
                result_seq.append("G")
            elif i=="A":
                result_seq.append("U")
            elif i=="U":
                result_seq.append("A")
            else:
                print("wrong RNA sequence")
            
        else:
            # DNA
            if i=="G":
                result_seq.append("C")
            elif i=="C":
                result_seq.append("G")
            elif i=="A":
                result_seq.append("T")
            elif i=="T":
                result_seq.append("A")
            else:
                print("wrong RNA sequence")
    return result_seq



def judge(piRNA,mRNA_sequence,pos):
    mRNA_start = int(pos.split('-')[0])-1
    mRNA_end = int(pos.split('-')[1])-1
    rem_seq=(mRNA_start,mRNA_end+1)

    pirna_com=complement(piRNA)
    target_score=[]

    if mRNA_start < 20:
        if mRNA_end+21 > len(mRNA_sequence):
            gu=gu_fun(0,len(mRNA_sequence)-21,mRNA_sequence,pirna_com,target_score,rem_seq,mRNA_start,mRNA_end)
        else:
            gu=gu_fun(0,mRNA_end,mRNA_sequence,pirna_com,target_score,rem_seq,mRNA_start,mRNA_end)
    elif mRNA_end+21 > len(mRNA_sequence):
        gu=gu_fun(mRNA_start-20,len(mRNA_sequence)-21,mRNA_sequence,pirna_com,target_score,rem_seq,mRNA_start,mRNA_end)
    else:
        gu=gu_fun(mRNA_start-20,mRNA_end,mRNA_sequence,pirna_com,target_score,rem_seq,mRNA_start,mRNA_end)
    return gu

