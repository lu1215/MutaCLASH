import pandas as pd
import time
import argparse
from collections import Counter
import re

parser = argparse.ArgumentParser()
parser.add_argument("--inputname", help="type input filename", type=str)

args = parser.parse_args()
inputname = args.inputname

data = pd.read_csv(inputname)
print(len(data))
start_m = [int(data['rem_tran_target_pos'][i].split('-')[0])+data['pos_m'][i]-1 for i in range(len(data))]
data['pos_m'] = start_m
del start_m

# target
D_dict = []
P_dict = []
M_dict = []
hybird_dict = []

time.sleep(1)
for i in range(len(data)):
    print(i, end='\r')
    read_count = data['read_count'][i]
    name = data['transcript_name'][i]
    D_list = []
    P_list = []
    M_list = []

    md = data['MD_m'][i][5::].replace('^', '')

    target_len = 0
    curr = int(data['pos_m'][i])
    for m in re.findall(r'[0-9]+|[A-Z]+',md):
        if m.isalpha():
            target_len += len(m)
        else:
            target_len += int(m)

    hybrid_list = '{}-{}'.format(curr, curr-1+target_len)

    curr = int(data['pos_m'][i])
    MD = data['MD_m'][i][5::]
    tmp = 0
    m = 0
    while(m<len(MD)):
        n = MD[m]
        if n.isalpha():
            if tmp != 0:
                P_list.append('{}-{}'.format(curr, curr+tmp-1))
            curr += tmp
            M_list.append(curr)
            curr += 1
            m += 1
            tmp = 0
        elif n == '^': 
            num_d = 0
            for d in range(m+1, len(MD)):
                if str(MD[d]).isalpha():
                    num_d += 1
                else:break

            if tmp != 0:
                P_list.append('{}-{}'.format(curr, curr+tmp-1))
            curr += tmp
            for t in range(num_d):
                D_list.append(curr+t)

            curr += num_d
            m += num_d+1
            tmp = 0
        else:
            tmp = 10*tmp + int(n)
            m += 1
    if tmp != 0:
        P_list.append('{}-{}'.format(curr, curr+tmp-1))
    D_dict.append(D_list)
    #P_dict.append(P_list)
    M_dict.append(M_list)
    hybird_dict.append(hybrid_list)
data['D'] = D_dict
data['M'] = M_dict
#data['target pos'] = hybird_dict
#data['A'] = [data['D'][i]+data['M'][i] for i in range(len(data))]
del D_list, P_list, M_list, hybrid_list
del P_dict, D_dict, M_dict, hybird_dict

#regulator
#     print(len(data))
#     D_dict = []
#     M_dict = []
#     time.sleep(1)
#     for i in trange(len(data)):
#         D_list = []
#         M_list = []
#         n = data['MD_r'][i][5:]
#         init_pos = data['pos_r'][i] - 1
#         tmp_pos = 0
#         m = 0
#         for j in n:

#             if j == '^':
#                 m += 1 + tmp_pos
#                 tmp_pos = 0
#                 D_list.append(init_pos+m)

#             elif j.isalpha():
#                 m += 1 + tmp_pos
#                 tmp_pos = 0
#                 M_list.append(init_pos+m)
#             else:
#                 tmp_pos = 10*tmp_pos + int(j)   

#     # print(D_list)
#     # print(M_list)
#         D_dict.append(D_list)
#         M_dict.append(M_list)
#     data['D_r'] = D_dict
#     data['M_r'] = M_dict

# map read count
data['count'] = 1
counter = Counter(data['hybrid0'])
readcount_nor = [data['read_count'][i]/counter[data['hybrid0'][i]] for i in range(len(data))]
count_nor = [data['count'][i]/counter[data['hybrid0'][i]] for i in range(len(data))]
data['nor_readcount'] = readcount_nor
data['nor_count'] = count_nor
del counter, readcount_nor, count_nor
outputname = inputname.split('/')[-1].replace('.csv', '')
data.to_csv('tmp/{}_detail.csv'.format(outputname), index=False)
