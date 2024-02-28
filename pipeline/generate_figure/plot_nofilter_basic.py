from plot_import import *

## mutation A、T、C、G
# 無突變、hybrid 上 target片段 ATCG比例
tmp = data[data['D'].isin(['[]'])]
tmp = tmp[tmp['M'].isin(['[]'])]
tmp = tmp.reset_index(drop=False)
a_num = 0
t_num = 0
c_num = 0
g_num = 0

for i in trange(len(tmp)):
    rc = tmp['nor_readcount'][i]
    region = tmp['remain_pos'][i].split('-')
    seq = tmp['hybrid_seq'][i][int(region[0])-1: int(region[1])]
    for m in seq:
        if m == 'A':
            a_num += rc
        elif m == 'T':
            t_num += rc
        elif m == 'G':
            g_num += rc
        elif m == 'C':
            c_num += rc
total_num = a_num+t_num+c_num+g_num
print('Without Mutation:')
print('a: {}, t: {}, c: {}, g: {}'.format(round(a_num/total_num, 2), round(t_num/total_num, 2), round(c_num/total_num, 2), round(g_num/total_num, 2)))

## deletion
datax = pd.merge(data, mrna_275, left_on='transcript_name', right_on='Gene name', how='inner')
mut = 'D'
tmp = datax[~datax[mut].isin(['[]'])]
tmp = tmp.reset_index(drop=False)
a_num = 0
t_num = 0
c_num = 0
g_num = 0

for i in trange(len(tmp)):
    rc = tmp['nor_readcount'][i]
    for m in eval(tmp[mut][i]):
        seq = 'sequence' # mRNA序列
        
        try:
            if tmp[seq][i][int(m)-1] == 'A': 
                a_num += rc
            elif tmp[seq][i][int(m)-1] == 'T':
                t_num += rc
            elif tmp[seq][i][int(m)-1] == 'G':
                g_num += rc
            elif tmp[seq][i][int(m)-1] == 'C':
                c_num += rc
        except:
            continue
            #print(int(m)+1)
            
total_num = a_num+t_num+c_num+g_num
print("\nDeletion:")
print('a: {}, t: {}, c: {}, g: {}'.format(round(a_num/total_num, 2), round(t_num/total_num, 2), round(c_num/total_num, 2), round(g_num/total_num, 2)))

## mismatch
datax = pd.merge(data, mrna_275, left_on='transcript_name', right_on='Gene name', how='inner')
mut = 'M'
tmp = datax[~datax[mut].isin(['[]'])]
tmp = tmp.reset_index(drop=False)
a_num = 0
t_num = 0
c_num = 0
g_num = 0

at_num = 0
ac_num = 0
ag_num = 0

ta_num = 0
tc_num = 0
tg_num = 0

ca_num = 0
ct_num = 0
cg_num = 0

ga_num = 0
gt_num = 0
gc_num = 0

for i in trange(len(tmp)):
    rc = tmp['nor_readcount'][i]
    for m in eval(tmp[mut][i]):
        seq = 'sequence' # sequence
        h = (int(m) - int(tmp['rem_tran_target_pos'][i].split('-')[0])) + int(tmp['remain_pos'][i].split('-')[0])
        for d in eval(tmp['D'][i]):
            if int(d) < int(m):
                h -= 1
        
        try:  
            hyb_seq = tmp['hybrid_seq'][i][h-1]
            s = tmp[seq][i][int(m)-1]
         
            if s == 'A':
                a_num += rc
                if hyb_seq == 'T':
                    at_num += rc
                elif hyb_seq == 'C':
                    ac_num += rc
                elif hyb_seq == 'G':
                    ag_num += rc
                
            elif s == 'T':
                t_num += rc
                if hyb_seq == 'A':
                    ta_num += rc
                elif hyb_seq == 'C':
                    tc_num += rc
                elif hyb_seq == 'G':
                    tg_num += rc
            elif s == 'G':
                g_num += rc
                if hyb_seq == 'A':
                    ga_num += rc
                elif hyb_seq == 'T':
                    gt_num += rc
                elif hyb_seq == 'C':   
                    gc_num += rc
                
            elif s == 'C':
                c_num += rc
                if hyb_seq == 'A':
                    ca_num += rc
                elif hyb_seq == 'T':
                    ct_num += rc
                elif hyb_seq == 'G':
                    cg_num += rc
                
        except:
            continue
total_num = a_num+t_num+c_num+g_num
a_total = at_num+ac_num+ag_num
t_total = ta_num+tc_num+tg_num
g_total = ga_num+gt_num+gc_num
c_total = ca_num+ct_num+cg_num

print("\nSubstitution:")
print('a: {}, at: {}, ac: {}, ag: {}'.format(round(a_num/total_num, 2),
                                         round(at_num/a_total, 2), round(ac_num/a_total, 2), round(ag_num/a_total, 2)))

print('t: {}, ta: {}, tc: {}, tg: {}'.format(round(t_num/total_num, 2),
                                         round(ta_num/t_total, 2), round(tc_num/t_total, 2), round(tg_num/t_total, 2)))

print('c: {}, ca: {}, ct: {}, cg: {}'.format(round(c_num/total_num, 2),
                                         round(ca_num/c_total, 2), round(ct_num/c_total, 2), round(cg_num/c_total, 2)))

print('g: {}, ga: {}, gt: {}, gc: {}'.format(round(g_num/total_num, 2),
                                         round(ga_num/g_total, 2), round(gt_num/g_total, 2), round(gc_num/g_total, 2)))

del datax

## mutation存在方式
'''
# 同時存在 deletion、mismatch
t = data[~data['D'].astype(str).isin(['[]'])]
t = t[~t['M'].astype(str).isin(['[]'])]
t = t[['D', 'M', 'read_count']]
print(sum(t['read_count']))
# one、two connect、two split 、three up
one = 0
two_connect = 0
two_split = 0
three_up = 0
s = 0
s_list = []
for i in trange(len(data)):
    n = data['A'][i]
    rc = data['nor_readcount'][i]
    l = eval(str(n))
    if len(l) == 1:
        one += rc
        
        
    elif len(l) == 2:
        s_list.append(s)
        if int(l[0])+1 == int(l[1]):
            two_connect += rc
        else:
            two_split += rc
    elif len(l) >= 3:
        three_up += rc
        s_list.append(s)
    s += 1
print(one)
print(two_connect)
print(two_split)
print(three_up)

one = 0
two_connect = 0
two_split = 0
three_up = 0
for i in trange(len(data)):
    n = data['M'][i]
    l = eval(str(n))
    rc = data['nor_readcount'][i]
    if len(l) == 1:
        one += rc
    elif len(l) == 2:
        if int(l[0])+1 == int(l[1]):
            two_connect += rc
        else:
            two_split += rc
    elif len(l) >= 3:
        three_up += rc
print('======= mismatch =======')
print(one)
print(two_connect)
print(two_split)
print(three_up)

one = 0
two_connect = 0
two_split = 0
three_up = 0
for i in trange(len(data)):
    n = data['D'][i]
    l = eval(str(n))
    rc = data['nor_readcount'][i]
    if len(l) == 1:
        one += rc
    elif len(l) == 2:
        if int(l[0])+1 == int(l[1]):
            two_connect += rc
        else:
            two_split += rc
    elif len(l) >= 3:
        three_up += rc
print('======= deletion =======')
print(one)
print(two_connect)
print(two_split)
print(three_up)
'''
