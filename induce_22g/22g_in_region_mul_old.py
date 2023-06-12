import pandas as pd
import argparse
import time
import multiprocessing
from tqdm import trange

def induce_region(idx, rem_tran_target_pos, transcript_name):
    tmp_out = {}
    wt_rc_sum = 0
    region = find_region_center(rem_tran_target_pos, transcript_name)
    ## wt
    try:
        tmp = wt_g22_data[wt_g22_data['ref_id'].isin([transcript_name])] 
        for i2, row2 in tmp.iterrows():
            row2_region = row2['rem_tran_target_pos'].split('-')
            start_tmp = int(row2_region[0])
            end_tmp = int(row2_region[1])
            start_cal = int(float(region.split('-')[0]))
            end_cal = int(float(region.split('-')[1]))

            if start_cal <= end_tmp and end_cal >= start_tmp:  # 100 region中有 22g
                wt_rc_sum += evenly_rc
            else:continue
    except:
        wt_rc_sum = 0
    ## mut
    #try:
    #    tmp = mut_g22_data[mut_g22_data['ref_id'].isin([transcript_name])] 
    #    for i2, row2 in tmp.iterrows():
    #        row2_region = rem_tran_target_pos.split('-')
    #        start_tmp = int(row2_region[0])
    #        end_tmp = int(row2_region[1])
    #        start_cal = int(float(region.split('-')[0]))
    #        end_cal = int(float(region.split('-')[1]))
    #        if start_cal <= end_tmp and end_cal >= start_tmp:  # 100 region中有 22g
    #            mut_rc_sum += evenly_rc
    #        else:continue
    #except:pass

    #tmp[idx] = {wt_rc_sum, mut_rc_sum}
    tmp_out[idx] = [idx, wt_rc_sum]
    return tmp_out

def induce_site(idx, rem_tran_target_pos, transcript_name):
    tmp_out = {}
    wt_rc_sum = 0
    mut_rc_sum = 0
    region = find_site_center(rem_tran_target_pos, transcript_name)
    ## wt
    try:
        tmp = wt_g22_data[wt_g22_data['ref_id'].isin([transcript_name])] 
        
        for i2, row2 in tmp.iterrows():
            row2_region = row2['rem_tran_target_pos'].split('-')
            start_tmp = int(row2_region[0])
            end_tmp = int(row2_region[1])
            start_cal = int(float(region.split('-')[0]))
            end_cal = int(float(region.split('-')[1]))
            evenly_rc = row2['evenly_rc']
            if start_cal <= end_tmp and end_cal >= start_tmp:  # 100 region中有 22g
                     
                wt_rc_sum += evenly_rc
            else:continue
    except:
        wt_rc_sum = 0
    # mut
    try:
        tmp = mut_g22_data[mut_g22_data['ref_id'].isin([transcript_name])] 
        for i2, row2 in tmp.iterrows():
            row2_region = row2['rem_tran_target_pos'].split('-')
            start_tmp = int(row2_region[0])
            end_tmp = int(row2_region[1])
            start_cal = int(float(region.split('-')[0]))
            end_cal = int(float(region.split('-')[1]))
            evenly_rc = row2['evenly_rc']
            if start_cal <= end_tmp and end_cal >= start_tmp:  # 100 region中有 22g
                mut_rc_sum += evenly_rc
            else:continue
    except:
        mut_rc_sum = 0

    tmp_out[idx] = [idx, wt_rc_sum, mut_rc_sum]
    return tmp_out 

def find_region_center(rem_tran_target_pos, transcript_name):
    target_pos = rem_tran_target_pos.split('-')
    center = (int(target_pos[0]) + int(target_pos[1]))/2 # identified region 找中心
    mrna_len = mrna_275_dict[transcript_name]
    if center%1 != 0:
        if center - 1 < ex_len:
            region = '{}-{}'.format(1, int(center)+ex_len+1)
        elif mrna_len - center < ex_len:
            region = '{}-{}'.format(int(center)-ex_len, mrna_len)
        else:
            region = '{}-{}'.format(int(center)-ex_len, int(center)+ex_len+1)
    else:
        if center - 1 < ex_len:
            region = '{}-{}'.format(1, center+ex_len)
        elif mrna_len - center < ex_len:
            region = '{}-{}'.format(center-ex_len, mrna_len)
        else:
            region = '{}-{}'.format(center-ex_len, center+ex_len)

    return region

def find_site_center(rem_tran_target_pos, transcript_name):
    target_pos = rem_tran_target_pos.split('-')
    center = (int(target_pos[0]) + int(target_pos[1]))/2 # site 找中心
    mrna_len = mrna_275_dict[transcript_name]
    if center%1 != 0:
        if center - 1 < ex_len:
            region = '{}-{}'.format(1, int(center)+ex_len+1)
        elif mrna_len - center < ex_len:
            region = '{}-{}'.format(int(center)-ex_len, mrna_len)
        else:
            region = '{}-{}'.format(int(center)-ex_len, int(center)+ex_len+1)
    else:
        if center - 1 < ex_len:
            region = '{}-{}'.format(1, center+ex_len)
        elif mrna_len - center < ex_len:
            region = '{}-{}'.format(center-ex_len, mrna_len)
        else:
            region = '{}-{}'.format(center-ex_len, center+ex_len)

    return region
parser = argparse.ArgumentParser()
parser.add_argument("--ex", help="type extend length", type=int)
parser.add_argument("--inputname", help="type CLASH result data", type=str)
parser.add_argument("--region_type", help="type region type", type=str)

args = parser.parse_args()
extend_len = args.ex
clash_result = args.inputname
region_type = args.region_type

print('==== data: {} ===='.format(clash_result))
# mRNA
mrna_275 = pd.read_csv('../data/reference/mRNA_WS275_WITHregion_v3.csv') 
mrna_275 = mrna_275[['Gene name', 'sequence']]
mrna_275_dict = {}
for i in range(len(mrna_275)):
    mrna_275_dict.update({mrna_275['Gene name'][i]: len(mrna_275['sequence'][i])})
del mrna_275
piRNA = pd.read_csv('../data/reference/piRNA_WS275.csv') #改rna檔案
piRNA_dict = {}
for i in range(len(piRNA)):
    piRNA_dict.update({piRNA['regulator_name'][i]: piRNA['raw_regulator_seq'][i]})

# 22G data
wt_g22_data = pd.read_csv('22g_data/WT_WAGOIP_EGL17M0_usemiRNA_norm_1step.csv')
mut_g22_data = pd.read_csv('22g_data/PRG_MUT_WAGOIP_EGL17M0_usemiRNA_norm_1step.csv')

# mRNA abundance
#m_abu = pd.read_csv('')
# clash data
clash_data = pd.read_csv(clash_result)
clash_data['raw_regulator_seq'] = [piRNA_dict[n.split('_')[0].replace(' ', '')] for n in clash_data['regulator_name']]
clash_data['idx'] = [i for i in range(len(clash_data))]
del piRNA_dict, piRNA

# code
ex_len = extend_len
wt_rc_list = []
mut_rc_list = []
data = clash_data.copy()
print(len(data))

if region_type == 'region':
    print(region_type)
    #data = data[['idx','transcript_name', 'rem_tran_target_pos']]
    
    cores = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(processes=cores-2) 
    b = {}
    d = [] 
    for i in trange(len(data)):

        wt_rc_sum = 0
        mut_rc_sum = 0
        transcript_name = data['transcript_name'][i]
        rem_tran_target_pos = data['rem_tran_target_pos'][i]
        idx = data['idx'][i]
        dicti = pool.apply_async(induce_region, (idx, rem_tran_target_pos, transcript_name))
        
        d.append(dicti)
    pool.close()
    pool.join()
    print('===== collecting result =====')
    print(len(d))
    for k in range(len(d)):
        items=d[k].get()
        b.update(items)        
    data_22g = pd.DataFrmae.from_dict(b, orient='index', columns=['idx', '22G_rc_WT', '22G_rc_MUT'])
    data = pd.merge(data, data_22g, on='idx', how='inner')
    print(len(data))
    data.to_csv('add_22g_info/22g_'+region_type+'_' + clash_result.split('/')[-1].split('.')[0] + '.csv', index=False)
    
if region_type == 'site':
    start = time.time()
    print(region_type)
    #data = data[['idx','transcript_name', 'rem_tran_target_pos', 'pirscan_target_endpos', 'raw_regulator_seq']]
    #data['pirscan_target_pos'] = ['{}-{}'.format(str(int(data['pirscan_target_endpos'][i])-len(data['raw_regulator_seq'][i])+1),
    #                                             str(int(data['pirscan_target_endpos'][i]))) for i in range(len(data))]

    cores = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(processes=cores-2) 
    b = {}
    d = [] 
    for i in trange(len(data)):

        transcript_name = data['transcript_name'][i]
        rem_tran_target_pos = data['pirscan_target_pos'][i]
        idx = data['idx'][i]
        dicti = pool.apply_async(induce_site, (idx, rem_tran_target_pos, transcript_name))
        d.append(dicti)
    pool.close()
    pool.join()
    print('===== collecting result =====')
    print(len(d))
    for k in range(len(d)):
        items=d[k].get()
        b.update(items)        
    data_22g = pd.DataFrame.from_dict(b, orient='index', columns=['idx', '22G_rc_WT', '22G_rc_MUT'])
    data = pd.merge(data, data_22g, on='idx', how='inner')
    data = data.drop(columns=['pirscan_target_pos'])
    print(len(data))
    stop = time.time()
    print((stop-start)/60)
    #print(data.head())
    data.to_csv('add_22g_info/22g_'+region_type+'_'+ str(ex_len) + '_' + clash_result.split('/')[-1].split('.')[0] + '.csv', index=False)
if region_type == 'abu':
    print(region_type)

