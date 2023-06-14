import pandas as pd
from ast import literal_eval
from tqdm import trange

data = '/media/disk1/ryan/tmp/WT1.csv'

df = pd.read_csv(data, usecols=['transcript_name','rem_tran_target_pos','M','D','count','read_count','nor_count','nor_readcount'],
                 converters={'M':literal_eval, 'D':literal_eval})

pos = df['rem_tran_target_pos'].str.split('-', expand=True).astype(int)
df['init_pos'] = pos[0]
df['end_pos'] = pos[1]
del df['rem_tran_target_pos']

D_lst, M_lst = [],[]

for i in trange(len(df)):
    if df.loc[i,'D']:
        D_lst.append(True)
    else:
        D_lst.append(False)
    if df.loc[i,'M']:
        M_lst.append(True)
    else:
        M_lst.append(False)

df['D_check'] = D_lst
df['M_check'] = M_lst

df_new = pd.read_csv('/media/disk1/ryan/tmp/WT1_new.csv',
                 converters={'hybrid_count':literal_eval,
                             'mismatch_count':literal_eval,
                             'deletion_count':literal_eval,
                             'hybrid_readcount':literal_eval,
                             'mismatch_readcount':literal_eval,
                             'deletion_readcount':literal_eval})


df_A = df[df['D_check'] | df['M_check']]
df_A['D_site1'] = [False]*len(df_A)
df_A['D_rc1'] = [False]*len(df_A)
df_A['M_site1'] = [False]*len(df_A)
df_A['M_rc1'] = [False]*len(df_A)

for i in trange(len(df_new)):
    mRNA = df_new.loc[i,'mRNA_name']
    df_tmp = df_A[ df_A['transcript_name']==mRNA ]
    
    for k,v in df_new.loc[i,'deletion_count'].items():
        if v<=1:
            lst = df_tmp[ df_tmp['D'].apply(lambda x: any([(v==k or v==str(k)) for v in x])) ].index.to_list()
            df_A.loc[df_A.index.isin(lst), 'D_site1'] = True

    for k,v in df_new.loc[i,'deletion_readcount'].items():
        if v<=1:
            lst = df_tmp[ df_tmp['D'].apply(lambda x: any([(v==k or v==str(k)) for v in x])) ].index.to_list()
            df_A.loc[df_A.index.isin(lst), 'D_rc1'] = True
    
    for k,v in df_new.loc[i,'mismatch_count'].items():
        if v<=1:
            lst = df_tmp[ df_tmp['M'].apply(lambda x: any([(v==k or v==str(k)) for v in x])) ].index.to_list()
            df_A.loc[df_A.index.isin(lst), 'M_site1'] = True
    
    for k,v in df_new.loc[i,'mismatch_readcount'].items():
        if v<=1:
            lst = df_tmp[ df_tmp['M'].apply(lambda x: any([(v==k or v==str(k)) for v in x])) ].index.to_list()
            df_A.loc[df_A.index.isin(lst), 'M_rc1'] = True


df_A.to_csv('/media/disk1/ryan/tmp/WT1_A.csv', index=False)
