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


def pvalue(f1,f2,outf):
    df1=pd.read_csv("up_"+f1,sep="\t",header=None)
    df2=pd.read_csv("up_"+f2,sep="\t",header=None)

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


    plt.figure()
    # plt.hist(data1,cumulative=True, density=1, bins=50,histtype='step',label=f1)
    # plt.hist(data2,cumulative=True, density=1, bins=50,histtype='step',label=f2)
    plt.hist(data1,cumulative=True, density=1, bins=len(data1),histtype='step',label=f1)
    plt.hist(data2,cumulative=True, density=1, bins=len(data2),histtype='step',label=f2)
    plt.legend(loc='upper left')
    plt.savefig('plot/hist_{}_{}.png'.format(prefix_f1,prefix_f2))


    outf.write(prefix_f1+"  "+prefix_f2+"\n")
    outf.write("-"*20+"\n")
    outf.write('%-20s %15.2f %15.2f %15.2f\n' %("KS_test",KStest[0],KStest[1],KStest[2]))
    outf.write('%-20s %15.2f %15.2f %15.2f\n' %("T_test",Ttest[0],Ttest[1],Ttest[2]))
    outf.write('%-20s %15.2f %15.2f %15.2f\n' %("U_test",Utest[0],Utest[1],Utest[2]))
    outf.write("\n\n")

def countUp(df1):
    groups = df1.groupby(df1[0])
    out = [[],[]]
    for name,group in groups:
        out[0].append(name)
        out[1].append(len(group))

    return out



def boxPlot(f1,f2,f3,jobID):

    df1=pd.read_csv("up_"+f1,sep="\t",header=None)
    df2=pd.read_csv("up_"+f2,sep="\t",header=None)
    df3=pd.read_csv("up_"+f3,sep="\t",header=None)
    res = pd.concat([df1,df2,df3],axis=1)
    res.columns=[f1,f2,f3]
    print("Filename:",jobID)
    print("mean:")
    print(res.mean())

    res.plot(kind='box', notch=False, grid=True)
    plt.title(jobID)
    plt.savefig('plot/boxPlot.png')


    c1 = countUp(df1)
    c2 = countUp(df2)
    c3 = countUp(df3)
    plt.figure()
    plt.scatter(c1[0],c1[1],s=5,alpha=0.5,c='r' ,label=f1)
    plt.scatter(c2[0],c2[1],s=5,alpha=0.5,c='y' ,label=f2)
    plt.scatter(c3[0],c3[1],s=5,alpha=0.5,c='b' ,label=f3)
    plt.title(jobID)
    plt.legend(loc='upper right')
    plt.savefig('plot/all_scatter.png')

    plt.figure()
    plt.subplot(2,2,1)
    plt.title(jobID)
    plt.ylabel(f1)
    plt.scatter(c1[0],c1[1],s=5,alpha=0.5,c='r')
    plt.subplot(222)
    plt.ylabel(f2)
    plt.scatter(c2[0],c2[1],s=5,alpha=0.5,c='y')
    plt.subplot(223)
    plt.ylabel(f3)
    plt.scatter(c3[0],c3[1],s=5,alpha=0.5,c='b')

    plt.savefig('plot/sep_scatter.png')

    plt.figure()
    plt.hist(df1.values,cumulative=True, density=1, bins=df1.shape[0],histtype='step',label=f1)
    plt.hist(df2.values,cumulative=True, density=1, bins=df2.shape[0],histtype='step',label=f2)
    plt.hist(df3.values,cumulative=True, density=1, bins=df3.shape[0],histtype='step',label=f3)
    plt.legend(loc='upper left')
    plt.title(jobID)
    
    plt.savefig('plot/hist_all.png')

def fileContent(f1):
    if os.path.getsize(f1):
        df1=pd.read_csv(f1,sep="\t",header=None)
        df1=df1.drop_duplicates(subset=None, keep='first', inplace=False)
        list1=df1.values.tolist()
    else:
        list1=[]
    return list1

def draw3(f1,f2,f3,cArray,jobID):
    prefix_f1=os.path.splitext(f1)[0]
    prefix_f2=os.path.splitext(f2)[0]
    prefix_f3=os.path.splitext(f3)[0]

    list1=fileContent(f1)
    list2=fileContent(f2)
    list3=fileContent(f3)

    jlist1=["/".join(i) for i in list1]
    jlist2=["/".join(i) for i in list2]
    jlist3=["/".join(i) for i in list3]
    print("set3",len(set(jlist1+jlist2+jlist3)))

    fig=plt.figure()
    venn3(subsets=[set(jlist1),set(jlist2),set(jlist3)],set_labels=[prefix_f1,prefix_f2,prefix_f3],set_colors=cArray)
    fig.suptitle(jobID)
    fig.savefig('plot/{}_{}_{}.png'.format(prefix_f1,prefix_f2,prefix_f3))
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


def draw2(f1,f2,cArray,jobID):
    prefix_f1=os.path.splitext(f1)[0]
    prefix_f2=os.path.splitext(f2)[0]

    list1=fileContent(f1)
    list2=fileContent(f2)


    jlist1=["".join(i) for i in list1]
    jlist2=["".join(i) for i in list2]

    fig=plt.figure()
    venn2(subsets=[set(jlist1),set(jlist2)],set_labels=[prefix_f1,prefix_f2],set_colors=cArray)
    fig.suptitle(jobID)
    fig.savefig('plot/{}_{}.png'.format(prefix_f1,prefix_f2))



def select3(way,jobID):
    df = pd.read_csv("/home/bba753951/Django/master_project/media/uploadfile/{}/{}/step6.csv".format(jobID,way),usecols=["hybrid_seq","regulator_name","transcript_name"])
    df.to_csv(way+".tab",index=0,header=0,sep="\t")


    df1 = pd.read_csv("/home/bba753951/Django/master_project/media/uploadfile/{}/{}/step6.csv".format(jobID,way),usecols=["RNAup_score"])
    df1.to_csv("up_"+way+".tab",index=0,header=0,sep="\t")
    

def main(number):

    dir_all=["pir","clan","hyb"]
    f1="pir.tab"
    f2="clan.tab"
    f3="hyb.tab"
    cArray=["red","#5599FF","green"]
    jobID=sys.argv[2]

    if number=="draw":
        draw3(f1,f2,f3,cArray,jobID)

        draw2(f1,f2,[cArray[0],cArray[1]],jobID)
        draw2(f1,f3,[cArray[0],cArray[2]],jobID)
        draw2(f3,f2,[cArray[2],cArray[1]],jobID)


    elif number=="drawall":
        drawAll(f1,f3,f2)
    elif number=="select":
        for i in dir_all:
            select3(i,jobID)
    elif number=="boxplot":
        boxPlot(f1,f2,f3,jobID)
    elif number=="pvalue":
        outfile = open("pvalue.txt","w")
        # outfile.write('%-20s %10s %10s %10s\n' %("pipeline","KS_test","T_test","U_test"))
        outfile.write("======= File name:{} ========\n".format(sys.argv[2]))
        pvalue(f1,f2,outfile)
        pvalue(f3,f2,outfile)
        pvalue(f1,f3,outfile)

        outfile.close()
    else:
        print("not found")


if __name__ == "__main__":
    main(sys.argv[1])

