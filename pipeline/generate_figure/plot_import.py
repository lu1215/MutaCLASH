import argparse
import pandas as pd
import numpy as np
import os
import re
import sys
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
    parser.add_argument("--trans", help="type trans data", type=str)
    parser.add_argument("--norm_factor", help="type base dataname", type=str)
    parser.add_argument("--tool", help="type base dataname", type=str)
    parser.add_argument("--region", help="type abundance regions", type=str)
    parser.add_argument("--figure", help="type file extension of figures", type=str)
    args = parser.parse_args()
    return args

def add_two_mRNA_list(new, target, gene):
    if gene == 0:
        return new
    elif gene == 1: # CSR-1 target
        csr1 = target[target['CSR1IP.N2__N2']==True].reset_index(drop=True) 
    elif gene == 2: # WAGO-1 target
        csr1 = target[target['WAGO1IP__WAGO1Input']==True].reset_index(drop=True)
    elif gene == 8: # Germline target
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
basename = args.basename
inputname = args.inputname
trans = args.trans
nor_f = float(args.norm_factor)
tool = args.tool
abu_region = args.region
fig_type = args.figure

## load data
d_name = basename
data = pd.read_csv(inputname)
mrna_275 = pd.read_csv(trans) 
mrna_275 = mrna_275[['Gene name', 'sequence']]

## Gene List
title_map_gene = {'0':'all mRNAs','1':'CSR-1 target','2':'WAGO-1 target', '8':'Germline target'}
