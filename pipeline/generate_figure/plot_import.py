import argparse
import pandas as pd
import numpy as np
import os
import re
from collections import Counter
import seaborn as sns
from tqdm import tqdm, trange
import time
import math
import scipy.stats as stats
from numpy import log as ln
import scipy.stats as stats
import scipy
import gc
import matplotlib
#matplotlib.use('agg')
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")
from statannot import add_stat_annotation

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--basename", help="type base dataname", type=str)
    parser.add_argument("--inputname", help="type input filename", type=str)
    parser.add_argument("--norm_factor", help="type base dataname", type=str)
    parser.add_argument("--tool", help="type base dataname", type=str)
    args = parser.parse_args()
    return args

target = pd.read_excel('../../data/reference/add_two_HCLee.RNAseq.master.xlsx')
def add_two_mRNA_list(new,gene):
    if gene == 0:
        return new
    elif gene == 1:
        csr1 = target[target['CSR1IP.N2__N2']==True].reset_index(drop=True) 
    elif gene == 2:
        csr1 = target[target['WAGO1IP__WAGO1Input']==True].reset_index(drop=True)
    elif gene == 8:
        csr1 = target[target['ce.germline.genes.Ortiz.G3_2014.type'].notnull()].reset_index(drop=True)
        
    tmp1 = new[new['Gene ID'].isin(list(csr1['row_names']))].reset_index(drop=True)
    return tmp1

def find_alpha(l):
    l2 = [i for i in l if i != 0]
    return min(l2)

def KS_test(x,y):
    if x and y:
        less = stats.mstats.ks_2samp(x, y,alternative = 'less')[1]
        greater = stats.mstats.ks_2samp(x, y,alternative = 'greater')[1]
      #  if (x.all() == y.all()) or (sum(x) == 0 and sum(y) == 0): #樣本相同 或 兩個樣本皆為0
      #      two_sided = 1.0
      #  else:
        two_sided = stats.mstats.ks_2samp(x, y,alternative = 'two-sided')[1]
    else:
        less = greater = two_sided = 0
    #print(two_sided, less, greater)
    return [two_sided, less, greater]

def T_test(x,y):
    if x and y:
        d, two_sided = stats.ttest_ind(x, y, equal_var=False)
        if d < 0:
            greater = 1 - two_sided/2 #"greater" is the alternative that x has a larger mean than y
            less = two_sided/2        #"less" is the alternative that x has a larger mean than y
        elif d >0:
            greater = two_sided/2
            less = 1 - two_sided/2
        else:
            greater = less = two_sided/2
    else:
        less = greater = two_sided = 0
    #print(two_sided, greater, less)
    return [two_sided, greater, less]

def U_test(x,y):
    if x and y:
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

    #print(two_sided, greater,less )
    return [two_sided, greater, less]

args = get_args()
inputname = args.inputname
basename = args.basename
nor_f = float(args.norm_factor)
tool = args.tool

mrna_275 = pd.read_csv('../../data/reference/mRNA_WS275_WITHregion_v3.csv') 
mrna_275 = mrna_275[['Gene name', 'sequence', 'Gene ID']]
mrna_275['Gene ID'] = mrna_275['Gene ID'].apply(lambda x:x.split('=')[1])

## Load Data
data = pd.read_csv(inputname)
d_name = basename

## user defined
title_map_gene = {'0':'all mRNAs','1':'CSR-1 target','2':'WAGO-1 target', '8':'Germline target'}
if tool=='pirScan':
    # targeting_score, mir_score, RNAup_score
    score_type =  'targeting_score'
    top = 10
    two_third = 0
    one_third = -15
    bot = -30
elif tool=='miRanda':
    score_type =  'mir_score'
    top = max(data[score_type])
    two_third = 140
    one_third = 100
    bot = 60
