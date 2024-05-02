import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib_venn import venn3
from matplotlib_venn import venn2
import sys

import scipy.stats as stats
def KS_test(x,y):
    if x.any() and y.any():
        less = stats.mstats.ks_2samp(x, y,alternative = 'less')[1] 
        greater = stats.mstats.ks_2samp(x, y,alternative = 'greater')[1]
      #  if (x.all() == y.all()) or (sum(x) == 0 and sum(y) == 0): #same sample or both equal 0
      #      two_sided = 1.0
      #  else:
        two_sided = stats.mstats.ks_2samp(x, y,alternative = 'two-sided')[1]
    else:
        less = greater = two_sided = 0


    return [two_sided, less, greater]
def T_test(x,y):
    if x.any() and y.any():
        
        d, two_sided = stats.ttest_ind(x, y, equal_var=False)
        if d < 0:
            greater = 1 - two_sided/2 #"greater" is the alternative that x has a larger mean than y
            less = two_sided/2 #"less" is the alternative that x has a larger mean than y
        elif d >0:
            greater = two_sided/2
            less = 1 - two_sided/2
        else:
            greater = less = two_sided/2
    else:
        less = greater = two_sided = 0

    return [two_sided, greater, less]
def U_test(x,y):
    if x.any() and y.any():
        d, two_sided = stats.ranksums(x, y)
        if d < 0:
            greater = 1 - two_sided/2
            less = two_sided/2
        elif d >0:
            greater = two_sided/2
            less = 1 - two_sided/2
        else:
            greater = less = two_sided/2
    else:
        less = greater = two_sided = 0
        
    return [two_sided, greater, less]


def pvalue(f1,f2,outf,kind):
    df1=pd.read_csv(kind+f1,sep="\t",header=None)
    df2=pd.read_csv(kind+f2,sep="\t",header=None)

    prefix_f1=os.path.splitext(f1)[0]
    prefix_f2=os.path.splitext(f2)[0]
    # print(prefix_f1,prefix_f2)
    

    data1=df1.iloc[:,0].values
    data2=df2.iloc[:,0].values

    KStest = KS_test(data1,data2)
    Ttest = T_test(data1,data2)
    Utest = U_test(data1,data2)

    # print("-"*20)
    # print(KStest)
    # print(Ttest)
    # print(Utest)
    # print("-"*20)


    # plt.hist(data1,cumulative=True, density=1, bins=50,histtype='step',label=f1)
    # plt.hist(data2,cumulative=True, density=1, bins=50,histtype='step',label=f2)

    # plt.figure()
    # plt.hist(data1,cumulative=True, density=1, bins=len(data1),histtype='step',label=f1)
    # plt.hist(data2,cumulative=True, density=1, bins=len(data2),histtype='step',label=f2)
    # plt.legend(loc='upper left')
    # plt.savefig('plot/hist_{}_{}.png'.format(prefix_f1,prefix_f2))


    outf.write(prefix_f1+"  "+prefix_f2+"\n")
    outf.write("-"*20+"\n")
    outf.write('%-20s %20s %20s %20s\n' %("KS_test",KStest[0],KStest[1],KStest[2]))
    outf.write('%-20s %20s %20s %20s\n' %("T_test",Ttest[0],Ttest[1],Ttest[2]))
    outf.write('%-20s %20s %20s %20s\n' %("U_test",Utest[0],Utest[1],Utest[2]))
    outf.write("\n\n")

def countUp(df1):
    groups = df1.groupby(df1[0])
    out = [[],[]]
    for name,group in groups:
        out[0].append(name)
        out[1].append(len(group))

    return out



def boxPlot(f1,f2,f3,jobID,kind):

    df1=pd.read_csv(kind+f1,sep="\t",header=None)
    df2=pd.read_csv(kind+f2,sep="\t",header=None)
    df3=pd.read_csv(kind+f3,sep="\t",header=None)
    res = pd.concat([df1,df2,df3],axis=1)
    f1=f1.replace(".tab","")
    f2=f2.replace(".tab","")
    f3=f3.replace(".tab","")
    # res.columns=[f1.replace(".tab",""),f2.replace(".tab",""),f3.replace(".tab","")]
    res.columns=[f1,f2,f3]
    print("Filename:",jobID)
    print("kind:",kind)
    print("mean:")
    print(res.mean())

    # res.plot(kind='box', notch=False, grid=True,showmeans=True)
    # plt.title(jobID)
    # plt.savefig('plot/{}boxPlot.png'.format(kind))


    f=res.boxplot(return_type='dict',patch_artist=True, notch=False, grid=True,showmeans=True,meanline=False,showbox=True,showfliers=True,whis=[0,100])
    for box in f['boxes']:
        box.set( color='black', linewidth=1 )
        box.set(facecolor='pink', alpha=0.8) 

    for whisker in f['whiskers']:
        whisker.set(color='k', linewidth=1, linestyle='--')

    for cap in f['caps']:
        cap.set(color='gray', linewidth=2)
    for median in f['medians']:
        median.set(color='orange', linewidth=2)
    for flier in f['fliers']:
        flier.set(marker='o', color='y', alpha=0.7)
    plt.title("{}({})".format(jobID,kind))
    plt.ylabel("RNAup sore")
    plt.savefig('plot/{}_boxPlot.png'.format(kind))
##############################################
    c1 = countUp(df1)
    c2 = countUp(df2)
    c3 = countUp(df3)
    plt.figure()
    plt.scatter(c1[0],c1[1],s=5,alpha=0.4,c='r' ,label=f1)
    plt.scatter(c2[0],c2[1],s=5,alpha=0.4,c='y' ,label=f2)
    plt.scatter(c3[0],c3[1],s=5,alpha=0.4,c='b' ,label=f3)
    plt.title(jobID)
    plt.legend(loc='upper right')
    plt.savefig('plot/{}all_scatter.png'.format(kind))

    plt.figure()
    plt.subplot(2,2,1)
    plt.title(jobID)
    plt.ylabel(f1)
    plt.scatter(c1[0],c1[1],s=5,alpha=0.4,c='r')
    plt.subplot(222)
    plt.ylabel(f2)
    plt.scatter(c2[0],c2[1],s=5,alpha=0.4,c='y')
    plt.subplot(223)
    plt.ylabel(f3)
    plt.scatter(c3[0],c3[1],s=5,alpha=0.4,c='b')

    plt.savefig('plot/{}sep_scatter.png'.format(kind))

###########################################
    plt.figure()
    plt.hist(df1.values,cumulative=True, density=1, bins=df1.shape[0],histtype='step',label=f1)
    plt.hist(df2.values,cumulative=True, density=1, bins=df2.shape[0],histtype='step',label=f2)
    plt.hist(df3.values,cumulative=True, density=1, bins=df3.shape[0],histtype='step',label=f3)
    plt.legend(loc='upper left')
    plt.title(jobID)
    
    plt.savefig('plot/{}hist_all.png'.format(kind))

def draw3(f1,f2,f3,cArray,jobID,kind):
    prefix_f1=os.path.splitext(f1)[0]
    prefix_f2=os.path.splitext(f2)[0]
    prefix_f3=os.path.splitext(f3)[0]

    df1=pd.read_csv(kind+f1,sep="\t",header=None)
    df1=df1.drop_duplicates(subset=None, keep='first', inplace=False)
    list1=df1.values.tolist()

    df2=pd.read_csv(kind+f2,sep="\t",header=None)
    df2=df2.drop_duplicates(subset=None, keep='first', inplace=False)
    list2=df2.values.tolist()

    df3=pd.read_csv(kind+f3,sep="\t",header=None)
    df3=df3.drop_duplicates(subset=None, keep='first', inplace=False)
    list3=df3.values.tolist()



    jlist1=["".join(i) for i in list1]
    jlist2=["".join(i) for i in list2]
    jlist3=["".join(i) for i in list3]
    print("set3",len(set(jlist1+jlist2+jlist3)))

    fig=plt.figure()
    venn3(subsets=[set(jlist1),set(jlist2),set(jlist3)],set_labels=[prefix_f1,prefix_f2,prefix_f3],set_colors=cArray)
    fig.suptitle(jobID)
    fig.savefig('plot/{}{}_{}_{}.png'.format(kind,prefix_f1,prefix_f2,prefix_f3))
    print("done")

def drawAll(f1,f2,f3):
    prefix_f1=os.path.splitext(f1)[0]
    prefix_f2=os.path.splitext(f2)[0]
    prefix_f3=os.path.splitext(f3)[0]

    df1=pd.read_csv(f1,sep="\t",header=None)
    df1=df1.drop_duplicates(subset=None, keep='first', inplace=False)
    list1=df1.values.tolist()

    df2=pd.read_csv(f2,sep="\t",header=None)
    df2=df2.drop_duplicates(subset=None, keep='first', inplace=False)
    list2=df2.values.tolist()

    df3=pd.read_csv(f3,sep="\t",header=None)
    df3=df3.drop_duplicates(subset=None, keep='first', inplace=False)
    list3=df3.values.tolist()



    jlist1=["".join(i) for i in list1]
    jlist2=["".join(i) for i in list2]
    jlist3=["".join(i) for i in list3]

    plt.figure()
    plt.subplot(2,2,1)
    venn2(subsets=[set(jlist1),set(jlist2)],set_labels=[prefix_f1,prefix_f2],set_colors=["red","green"], alpha = 0.4)
    plt.subplot(2,2,2)
    venn2(subsets=[set(jlist2),set(jlist3)],set_labels=[prefix_f2,prefix_f3],set_colors=["green","#5599FF"], alpha = 0.4)
    plt.subplot(2,2,3)
    venn2(subsets=[set(jlist1),set(jlist3)],set_labels=[prefix_f1,prefix_f3],set_colors=["red","#5599FF"], alpha = 0.4)
    plt.subplot(2,2,4)
    venn3(subsets=[set(jlist1),set(jlist2),set(jlist3)],set_labels=[prefix_f1,prefix_f2,prefix_f3],set_colors=["red","green","#5599FF"], alpha = 0.4)
    plt.savefig('plot/together_{}_{}_{}.png'.format(prefix_f1,prefix_f2,prefix_f3))
    print("done")


def draw2(f1,f2,cArray,jobID,kind):
    prefix_f1=os.path.splitext(f1)[0]
    prefix_f2=os.path.splitext(f2)[0]


    df1=pd.read_csv(kind+f1,sep="\t",header=None)
    df1=df1.drop_duplicates(subset=None, keep='first', inplace=False)
    list1=df1.values.tolist()

    df2=pd.read_csv(kind+f2,sep="\t",header=None)
    df2=df2.drop_duplicates(subset=None, keep='first', inplace=False)
    list2=df2.values.tolist()


    jlist1=["".join(i) for i in list1]
    jlist2=["".join(i) for i in list2]

    fig=plt.figure()
    venn2(subsets=[set(jlist1),set(jlist2)],set_labels=[prefix_f1,prefix_f2],set_colors=cArray)
    fig.suptitle(jobID)
    fig.savefig('plot/{}{}_{}.png'.format(kind,prefix_f1,prefix_f2))

def pAll(data1,data2):
    KStest = KS_test(data1,data2)
    Ttest = T_test(data1,data2)
    Utest = U_test(data1,data2)
    # print("KS_test\t",KStest)
    # print("T_test\t",Ttest)
    # print("U_test\t",Utest)
    print('%-20s %20s %20s %20s\n' %("KS_test",KStest[0],KStest[1],KStest[2]))
    print('%-20s %20s %20s %20s\n' %("T_test",Ttest[0],Ttest[1],Ttest[2]))
    print('%-20s %20s %20s %20s\n' %("U_test",Utest[0],Utest[1],Utest[2]))

def affectBar(len_data,data_count,kind,jobID,filt,columns,kind2):
    plt.figure(figsize=(10,4))
    Y=np.array(len_data)
    X=np.arange(data_count)
    plt.bar(X,Y,width = 0.5, facecolor = '#009FCC', edgecolor= 'black')
    alll=len_data[0]
    for i, j in zip(X, Y):
        plt.text(i-0.2, j+alll*0.01, '%d'%j, color='red')
        plt.text(i-0.2, j/2, '%.2f %s'%(j/alll*100,"%"), color='white')
    plt.title("{}({})".format(jobID,kind))
    plt.xlabel("{}({} equal)".format(kind,filt))
    plt.ylabel("remain {} count".format(kind2))
    plt.xticks(range(data_count),columns)
    plt.savefig('plot/{}_bar_{}.png'.format(kind,kind2))

    print(kind2)
    for i in range(data_count):
        print("%s(%.2f"%(len_data[i],len_data[i]/alll*100))
        if i == data_count - 1:
            break

def ori_affectBar(len_data,ori_data,data_count,kind,jobID,filt,columns,kind2):
    plt.figure(figsize=(10,4))
    np_len=np.array(len_data)
    np_ori=np.array(ori_data)
    print("\n\nori_select",ori_data)
    print("hits_select",len_data)
    Y=np_len/np_ori*100
    X=np.arange(data_count)
    plt.bar(X,Y,width = 0.5, facecolor = '#009FCC', edgecolor= 'black')
    alll=Y[0]
    for i, j in zip(X, Y):
        plt.text(i-0.2, j+alll*0.01, '%.2f %s'%(j,"%"), color='red')
    plt.title("{}({})".format(jobID,kind))
    plt.xlabel("{}({} equal)".format(kind,filt))
    plt.ylabel("hits percent")
    plt.xticks(range(data_count),columns)
    plt.savefig('plot/hits_{}_bar_{}.png'.format(kind,kind2))

    print("------------hits persent of {}-----------".format(kind2))
    print(kind2)
    fp=open("hits_{}".format(kind),"w")
    fp.write(jobID[-2:]+",")
    for i in range(data_count):
        print("%.3f"%(Y[i]))
        if i < data_count-1:
            fp.write("%.3f,"%(Y[i]))
        else:
            fp.write("%.3f\n"%(Y[i]))

    fp.close()

        


#kind
# 0 for fold
# 1 for readcount
def Affect(way,jobID,filterArray,kind):
    df = pd.read_csv("/home/bba753951/Django/master_project/media/uploadfile/{}/{}/hyb_file_step5.csv".format(jobID,way),usecols=["hybrid_seq","regulator_name","transcript_name","RNAup_pos","RNAup_score","RNAfold_MFE","read_count"])
    step2df = pd.read_csv("/home/bba753951/Django/master_project/media/uploadfile/{}/hyb_file_step2.csv".format(jobID),usecols=["RNAfold_MFE","read_count"])
    # df=df.drop_duplicates(subset=None, keep='first', inplace=False)
    print("original step2 count",len(step2df))
    if kind=="read_count":
        filt="greater"
        df1=df.loc[:,["read_count","RNAup_score","hybrid_seq"]]
        fold0=df1["read_count"] >= filterArray[0]
        fold5=df1["read_count"] >= filterArray[1]
        fold10=df1["read_count"] >= filterArray[2]
        fold15=df1["read_count"] >= filterArray[3]
        fold20=df1["read_count"] >= filterArray[4]

        ori0=step2df["read_count"] >= filterArray[0]
        ori5=step2df["read_count"] >= filterArray[1]
        ori10=step2df["read_count"] >= filterArray[2]
        ori15=step2df["read_count"] >= filterArray[3]
        ori20=step2df["read_count"] >= filterArray[4]

    elif kind=="RNAfold":
        filt="less"
        df1=df.loc[:,["RNAfold_MFE","RNAup_score","hybrid_seq"]]
        fold0=df1["RNAfold_MFE"] <= filterArray[0]
        fold5=df1["RNAfold_MFE"] <= filterArray[1]
        fold10=df1["RNAfold_MFE"] <= filterArray[2]
        fold15=df1["RNAfold_MFE"] <= filterArray[3]
        fold20=df1["RNAfold_MFE"] <= filterArray[4]

        ori0=step2df["RNAfold_MFE"] <= filterArray[0]
        ori5=step2df["RNAfold_MFE"] <= filterArray[1]
        ori10=step2df["RNAfold_MFE"] <= filterArray[2]
        ori15=step2df["RNAfold_MFE"] <= filterArray[3]
        ori20=step2df["RNAfold_MFE"] <= filterArray[4]
    else:
        print("error kind")
        raise ValueError

    pdata_f15=df1[(fold15)]
    pdata_f10=df1[(fold10)]
    pdata_f20=df1[(fold20)]
    pdata_f5=df1[(fold5)]
    pdata_f0=df1[(fold0)]


    count_ori_f15=len(step2df[(ori15)])
    count_ori_f10=len(step2df[(ori10)])
    count_ori_f20=len(step2df[(ori20)])
    count_ori_f5=len(step2df[(ori5)])
    count_ori_f0=len(step2df[(ori0)])

    data_f15=pdata_f15.loc[:,"RNAup_score"]
    data_f10=pdata_f10.loc[:,"RNAup_score"]
    data_f20=pdata_f20.loc[:,"RNAup_score"]
    data_f5=pdata_f5.loc[:,"RNAup_score"]
    data_f0=pdata_f0.loc[:,"RNAup_score"]






    hyb_f0=pdata_f0[["hybrid_seq"]]
    hyb_f0=hyb_f0.drop_duplicates(subset=None, keep='first', inplace=False)
    hyb_f5=pdata_f5[["hybrid_seq"]]
    hyb_f5=hyb_f5.drop_duplicates(subset=None, keep='first', inplace=False)
    hyb_f10=pdata_f10[["hybrid_seq"]]
    hyb_f10=hyb_f10.drop_duplicates(subset=None, keep='first', inplace=False)
    hyb_f15=pdata_f15[["hybrid_seq"]]
    hyb_f15=hyb_f15.drop_duplicates(subset=None, keep='first', inplace=False)
    hyb_f20=pdata_f20[["hybrid_seq"]]
    hyb_f20=hyb_f20.drop_duplicates(subset=None, keep='first', inplace=False)

    all_hyb=[hyb_f0,hyb_f5,hyb_f10,hyb_f15,hyb_f20]
    len_hyb=list(map(len,all_hyb))
    hyb_count=len(all_hyb)

    all_data=[data_f0,data_f5,data_f10,data_f15,data_f20]
    len_data=list(map(len,all_data))
    data_count=len(all_data)
    columns=list(map(str,filterArray))

    len_ori=[]

    all_count_ori=[count_ori_f0,count_ori_f5,count_ori_f10,count_ori_f15,count_ori_f20]

    affectBar(len_data,data_count,kind,jobID,filt,columns,"pair")
    # remain read which can make pair
    affectBar(len_hyb,hyb_count,kind,jobID,filt,columns,"sequence")
    # remain original CLASH read
    affectBar(all_count_ori,hyb_count,kind,jobID,filt,columns,"\'CLASH read\'")

    # ori_affectBar(len_data,all_count_ori,data_count,kind,jobID,filt,columns,"pair")
    ori_affectBar(len_hyb,all_count_ori,hyb_count,kind,jobID,filt,columns,"sequence")

    plt.figure(figsize=(10,4))
    res = pd.concat(all_data,axis=1)
    res.columns=columns
    print("\n{} mean:".format(kind))
    print(res.mean())
    print("\n{} pvalue:".format(kind))
    for i in range(data_count):
        print("\n{} vs {}".format(columns[i],columns[i+1]))
        pAll(all_data[i],all_data[i+1])
        if i == data_count-2:
            break


    f=res.boxplot(return_type='dict',patch_artist=True, notch=False, grid=True,showmeans=True,meanline=False,showbox=True,showfliers=True,whis=[0,100])
    for box in f['boxes']:
        box.set( color='black', linewidth=1 )
        box.set(facecolor='pink', alpha=0.8) 

    for whisker in f['whiskers']:
        whisker.set(color='k', linewidth=1, linestyle='--')

    for cap in f['caps']:
        cap.set(color='gray', linewidth=2)
    for median in f['medians']:
        median.set(color='orange', linewidth=2)
    for flier in f['fliers']:
        flier.set(marker='o', color='y', alpha=0.7)
    plt.title("{}({})".format(jobID,kind))
    plt.xlabel("{}({} equal)".format(kind,filt))
    plt.ylabel("RNAup sore")
    plt.savefig('plot/{}_boxPlot.png'.format(kind))

    plt.figure()
    plt.hist(data_f0.values,cumulative=True, density=1, bins=data_f0.shape[0],histtype='step',label=columns[0])
    plt.hist(data_f5.values,cumulative=True, density=1, bins=data_f5.shape[0],histtype='step',label=columns[1])
    plt.hist(data_f10.values,cumulative=True, density=1, bins=data_f10.shape[0],histtype='step',label=columns[2])
    plt.hist(data_f15.values,cumulative=True, density=1, bins=data_f15.shape[0],histtype='step',label=columns[3])
    plt.hist(data_f20.values,cumulative=True, density=1, bins=data_f20.shape[0],histtype='step',label=columns[4])
    plt.legend(loc='upper left')
    plt.title("{}({})".format(jobID,kind))
    plt.savefig('plot/hist_affect_{}.png'.format(kind))

def select3(way,jobID):
    df = pd.read_csv("/home/bba753951/Django/master_project/media/uploadfile/{}/{}/hyb_file_step5.csv".format(jobID,way),usecols=["hybrid_seq","regulator_name","transcript_name","RNAup_pos","RNAup_score"])
    # df=df.drop_duplicates(subset=None, keep='first', inplace=False)

    df1=df.loc[:,["hybrid_seq","regulator_name","transcript_name"]]
    df1.to_csv(way+".tab",index=0,header=0,sep="\t")
    # print(df1.head(3))
    df1=df.loc[:,["regulator_name","transcript_name"]]
    df1.to_csv("pair_"+way+".tab",index=0,header=0,sep="\t")


    # df1 = pd.read_csv("/home/bba753951/Django/master_project/media/uploadfile/{}/{}/hyb_file_step5.csv".format(jobID,way),usecols=["RNAup_score"])
    df1=df.loc[:,["RNAup_score"]]
    # print(df1.head(3))
    df1.to_csv("up_"+way+".tab",index=0,header=0,sep="\t")

    df2 = pd.read_csv("/home/bba753951/Django/master_project/media/uploadfile/{}/{}/hyb_file_step5.csv".format(jobID,way),usecols=["hybrid_seq","RNAfold_MFE"])
    # print(df2.shape)
    df2=df2.drop_duplicates(subset=None, keep='first', inplace=False)
    # print(df2.shape)
    # print(df2.head(3))
    df3=df2.loc[:,["RNAfold_MFE"]]
    df3.to_csv("fold_"+way+".tab",index=0,header=0,sep="\t")
    

def main(number):

    dir_all=["pir","clan","hyb"]
    f1="pir.tab"
    f2="clan.tab"
    f3="hyb.tab"
    cArray=["red","#5599FF","green"]
    jobID=sys.argv[2]

    if number=="draw":
        draw3(f1,f2,f3,cArray,jobID,"")

        draw2(f1,f2,[cArray[0],cArray[1]],jobID,"")
        draw2(f1,f3,[cArray[0],cArray[2]],jobID,"")
        draw2(f3,f2,[cArray[2],cArray[1]],jobID,"")

        draw3(f1,f2,f3,cArray,jobID,"pair_")

        draw2(f1,f2,[cArray[0],cArray[1]],jobID,"pair_")
        draw2(f1,f3,[cArray[0],cArray[2]],jobID,"pair_")
        draw2(f3,f2,[cArray[2],cArray[1]],jobID,"pair_")

    elif number=="drawall":
        drawAll(f1,f3,f2)
    elif number=="affect":
        Affect("clan",jobID,[0,-5,-10,-15,-20],"RNAfold")
        Affect("clan",jobID,[1,2,3,4,5],"read_count")
    elif number=="select":
        for i in dir_all:
            select3(i,jobID)
    elif number=="boxplot":
        boxPlot(f1,f2,f3,jobID,"up_")
        # boxPlot(f1,f2,f3,jobID,"fold_")
    elif number=="pvalue":
        outfile = open("up_pvalue.txt","w")
        outfile.write("======= File name:{} ========\n".format(sys.argv[2]))
        pvalue(f1,f2,outfile,"up_")
        pvalue(f3,f2,outfile,"up_")
        pvalue(f1,f3,outfile,"up_")

        outfile.close()

        outfile = open("fold_pvalue.txt","w")
        outfile.write("======= File name:{} ========\n".format(sys.argv[2]))
        pvalue(f1,f2,outfile,"fold_")
        pvalue(f3,f2,outfile,"fold_")
        pvalue(f1,f3,outfile,"fold_")

        outfile.close()
    else:
        print("not found")


if __name__ == "__main__":
    main(sys.argv[1])

