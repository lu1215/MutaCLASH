import pandas as pd
import sys
import os
def splitPair(f1,f2):
    df=pd.read_csv(f1,sep="\t",header=None)
    df=df.drop_duplicates(subset=None, keep='first', inplace=False)
    groups = df.groupby(df[0])

    output=[]
    data = pd.DataFrame()
    for name,group in groups:
        if len(group) >1:
            groups2 = group.groupby(df[1])
            for name2,group2 in groups2:
                if len(group2) >1:
                    seq1=group2.iloc[0,5].split(";")
                    seq2=group2.iloc[1,5].split(";")
                    pair = [[a, b] for a in seq1  
                                    for b in seq2] 

                    # pair_p = [[i[0].split(":")[0],i[1].split(":")[0],i[0].split(":")[1],i[1].split(":")[1]] for i in pair ]

                    hyb_pos1="-".join(map(str,group2.iloc[0,2:4]))
                    hyb_pos2="-".join(map(str,group2.iloc[1,2:4]))
                    hyb_name=group2.iloc[0,0]
                    for i in pair:
                        nameID1,pos1=i[0].split(":")
                        nameID2,pos2=i[1].split(":")
                        if "tran" in nameID1:
                            ori_col=[hyb_name,nameID2,hyb_pos2,pos2,nameID1,hyb_pos1,pos1]
                        else:
                            ori_col=[hyb_name,nameID1,hyb_pos1,pos1,nameID2,hyb_pos2,pos2]
                        output.append(ori_col)


    data = data.append(output)
    data.columns = ["hyb_seq","regulator0","reg_hyb_target_pos","on_reg_pos","transcript0","remain_pos","rem_tran_target_pos"]
    print(data.head())
    data.to_csv(f2,sep=",",index=0,header=1)


def removePair_pos(f1,f2,n):
    n=int(n)
    df1=pd.read_csv(f1,index_col="hybrid_seq")
    print(df1.head())
    print("original:",len(df1))
    groups = df1.groupby(df1.index)
    count=0
    a=1
    
    for name,group in groups:
        if len(group) > (n/2-1) :
            count_set=set()
            list1=group[["regulator0","reg_hyb_target_pos","transcript0","remain_pos"]].values.tolist()
            for i in list1:
                count_set.update([i[0]+":"+i[1],i[2]+":"+i[3]])


            if len(count_set) > n:
                count+=len(group)
                df1=df1.drop(index=[name],axis=0)
                if a==1:
                    group.to_csv("bbb.csv",header=1,index=1)
                    a=0
        

    print("remain:",len(df1))
    print("remove:",count)

    # add column for match hyb and clan
    df1['on_reg_pos']=df1.apply(lambda x: "1-"+str(x['regulator_len']),axis=1)

    df1.to_csv(f2,header=1,index=1)


def removePair(f1,f2,n):
    n=int(n)
    df1=pd.read_csv(f1,index_col="hybrid_seq")
    print(df1.head())
    print("original:",len(df1))
    groups = df1.groupby(df1.index)
    count=0
    
    for name,group in groups:
        if len(group) > (n/2-1) :
            count_set=set()
            list1=group[["regulator0","transcript0"]].values.tolist()
            for i in list1:
                count_set.update(i)


            if len(count_set) > n:
                count+=len(group)
                df1=df1.drop(index=[name],axis=0)
        

    print("remain:",len(df1))
    print("remove:",count)

    # add column for match hyb and clan
    df1['on_reg_pos']=df1.apply(lambda x: "1-"+str(x['regulator_len']),axis=1)

    df1.to_csv(f2,header=1,index=1)


def main(way):
    if way == "splitPair":
        splitPair(sys.argv[2],sys.argv[3])
    elif way == "removePair":
        removePair(sys.argv[2],sys.argv[3],sys.argv[4])
    elif way == "removePair_pos":
        removePair_pos(sys.argv[2],sys.argv[3],sys.argv[4])
    else:
        print("No way = ",way)


if __name__ == "__main__":
    main(sys.argv[1])
