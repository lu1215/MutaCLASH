from plot_import import *

if tool=='pirScan':
    score_type =  'targeting_score'
elif tool=='miRanda':
    score_type =  'mir_score'
elif tool=='RNAup':
    score_type =  'RNAup_score'
    data = data[data[score_type] != 100]
    print("Filter Score != 100:{}".format(len(data)))
else:
    sys.exit(0)

if abu_region:
    top, two_third, one_third, bot = [float(x) for x in abu_region.split('/')]
else:
    one_third = np.quantile(list(data[score_type]), 0.333)
    two_third = np.quantile(list(data[score_type]), 0.666)
    top = max(data[score_type])
    bot = min(data[score_type])
    
print('high: {} < S <= {}'.format(str(two_third), str(top)))
print('mid: {} < S <= {}'.format(str(one_third), str(two_third)))
print('low: {} < S <= {}'.format(str(bot), str(one_third)))

# Gene List
title_map_gene = {'0':'all mRNAs','1':'CSR-1 target','2':'WAGO-1 target', '8':'Germline target'}
target = pd.read_excel('../../data/reference/add_two_HCLee.RNAseq.master.xlsx')
if 'Gene ID' not in data:
    print('no "Gene ID", using reference file: mRNA_WS275_IDtoName.csv')
    data_id = pd.read_csv('../../data/reference/mRNA_WS275_IDtoName.csv')
    data_id = data_id.rename(columns={'Gene name':'transcript_name'})
    data = pd.merge(data, data_id, on='transcript_name', how='inner')

# fold change
data = data.fillna(0)

alpha1 = find_alpha(data['alg1.r1.rpkm'])
alpha2 = find_alpha(data['alg1.r2.rpkm'])
alpha3 = find_alpha(data['alg1.r3.rpkm'])
alpha4 = find_alpha(data['wt.r1.rpkm'])
alpha5 = find_alpha(data['wt.r2.rpkm'])
alpha6 = find_alpha(data['wt.r3.rpkm'])
alpha = min([alpha1, alpha2, alpha3, alpha4, alpha5, alpha6])
print('alpha: {}'.format(alpha))

data['alg1.r1.rpkm_alpha'] = [data['alg1.r1.rpkm'][i]+alpha for i in range(len(data))]
data['alg1.r2.rpkm_alpha'] = [data['alg1.r2.rpkm'][i]+alpha for i in range(len(data))]
data['alg1.r3.rpkm_alpha'] = [data['alg1.r3.rpkm'][i]+alpha for i in range(len(data))]
data['wt.r1.rpkm_alpha'] = [data['wt.r1.rpkm'][i]+alpha for i in range(len(data))]
data['wt.r2.rpkm_alpha'] = [data['wt.r2.rpkm'][i]+alpha for i in range(len(data))]
data['wt.r3.rpkm_alpha'] = [data['wt.r3.rpkm'][i]+alpha for i in range(len(data))]

data['wt_avg_rpkm'] = [(data['wt.r1.rpkm'][i]+data['wt.r2.rpkm'][i]+data['wt.r3.rpkm'][i])/3 for i in range(len(data))]
data['alg1_avg_rpkm'] = [(data['alg1.r1.rpkm'][i]+data['alg1.r2.rpkm'][i]+data['alg1.r3.rpkm'][i])/3 for i in range(len(data))]

data['wt_avg_rpkm_alpha'] = [data['wt_avg_rpkm'][i]+alpha for i in range(len(data))]
data['alg1_avg_rpkm_alpha'] = [data['alg1_avg_rpkm'][i]+alpha for i in range(len(data))]
data['fold_change_avg'] = [math.log(data['alg1_avg_rpkm_alpha'][i]/data['wt_avg_rpkm_alpha'][i], 2) for i in range(len(data))]

data['fold_change1'] = [math.log(data['alg1.r1.rpkm_alpha'][i]/data['wt.r1.rpkm_alpha'][i], 2) for i in range(len(data))]
data['fold_change2'] = [math.log(data['alg1.r2.rpkm_alpha'][i]/data['wt.r2.rpkm_alpha'][i], 2) for i in range(len(data))]
data['fold_change3'] = [math.log(data['alg1.r3.rpkm_alpha'][i]/data['wt.r3.rpkm_alpha'][i], 2) for i in range(len(data))]
print('number: {}'.format(len(data)))

x1 = [n for n in data['alg1.r1.rpkm']]
y1 = [n for n in data['wt.r1.rpkm']]
no_0_list = []
for i,j in zip(x1, y1):
    if i != 0 and j != 0:
        no_0_list.append(math.log(i/j, 2))
    else:
        no_0_list.append("NULL")
data['fold_change1_without0'] = no_0_list
tmp_list = no_0_list.copy()
tmp_list.remove("NULL")
print('number of fold-change1 (without0): {}'.format(len(tmp_list)))

x1 = [n for n in data['alg1.r2.rpkm']]
y1 = [n for n in data['wt.r2.rpkm']]
no_0_list = []
for i,j in zip(x1, y1):
    if i != 0 and j != 0:
        no_0_list.append(math.log(i/j, 2))
    else:
        no_0_list.append("NULL")
data['fold_change2_without0'] = no_0_list
tmp_list = no_0_list.copy()
tmp_list.remove("NULL")
print('number of fold-change2 (without0): {}'.format(len(tmp_list)))

x1 = [n for n in data['alg1.r3.rpkm']]
y1 = [n for n in data['wt.r3.rpkm']]
no_0_list = []
for i,j in zip(x1, y1):
    if i != 0 and j != 0:
        no_0_list.append(math.log(i/j, 2))
    else:
        no_0_list.append("NULL")
data['fold_change3_without0'] = no_0_list
tmp_list = no_0_list.copy()
tmp_list.remove("NULL")
print('number of fold-change3 (without0): {}'.format(len(tmp_list)))

x1 = [n for n in data['alg1_avg_rpkm']]
y1 = [n for n in data['wt_avg_rpkm']]
no_0_list = []
for i,j in zip(x1, y1):
    if i != 0 and j != 0:
        no_0_list.append(math.log(i/j, 2))
    else:
        no_0_list.append("NULL")
data['fold_change_avg_without0'] = no_0_list
tmp_list = no_0_list.copy()
tmp_list.remove("NULL")
print('number of average fold-change (without0): {}\n'.format(len(tmp_list)))

# mRNA abundance fold change
dn_list = ['all data', 'high score', 'middle score', 'low score']

for group in [0]:
    for mut in ['D', 'M']:
        abu_num = 0
        for abu in ['fold_change_avg']:
            abu_num += 1
        
            if mut == 'D':
                mut_n = 'del'
            elif mut == 'M':
                mut_n = 'mis'
            elif mut == 'A':
                mut_n = 'mut'

            ## compare with all mutation (contorol group)
            print('======= compare with all mutation =======')
            # split group
            ana_data = add_two_mRNA_list(data, target, group)
            no_del_clash_result = ana_data[ana_data['A'].astype(str).isin(['[]'])]
            del_clash_result = ana_data[~ana_data[mut].astype(str).isin(['[]'])] 
            ## 無突變資料中將突變資料的 transcript去除 (看有沒有需要)
            no_del_clash_result = no_del_clash_result[~no_del_clash_result['transcript_name'].isin(list(del_clash_result['transcript_name']))]
            ####
            print(len(del_clash_result), len(no_del_clash_result))
            # three different score
            no_one_third_data = no_del_clash_result[no_del_clash_result[score_type]<=one_third]
            no_one_third_data = no_one_third_data[no_one_third_data[score_type]>bot]
            no_two_third_data = no_del_clash_result[no_del_clash_result[score_type]>one_third]
            no_two_third_data = no_two_third_data[no_two_third_data[score_type]<=two_third]
            no_three_third_data = no_del_clash_result[no_del_clash_result[score_type]>two_third]
            no_three_third_data = no_three_third_data[no_three_third_data[score_type]<=top]

            del_one_third_data = del_clash_result[del_clash_result[score_type]>bot]
            del_one_third_data = del_one_third_data[del_one_third_data[score_type]<=one_third]
            del_two_third_data = del_clash_result[del_clash_result[score_type]>one_third]
            del_two_third_data = del_two_third_data[del_two_third_data[score_type]<=two_third]
            del_three_third_data = del_clash_result[del_clash_result[score_type]>two_third]
            del_three_third_data = del_three_third_data[del_three_third_data[score_type]<=top]

            d1_list = [del_clash_result, del_three_third_data, del_two_third_data, del_one_third_data]
            d2_list = [no_del_clash_result, no_three_third_data, no_two_third_data, no_one_third_data]
            d_n = 0
#             print('no {}: '.format('mut'), len(no_three_third_data), len(no_two_third_data), len(no_one_third_data))
#             print('{}: '.format(mut), len(del_three_third_data), len(del_two_third_data), len(del_one_third_data))

            for d1, d2 in zip(d1_list, d2_list):

                print(d_n, title_map_gene[str(group)], dn_list[d_n], mut_n, abu, 'with mutation')
                # plot
                x = d1[abu]
                y = d2[abu]
                tmp1 = pd.Series(list(x))
                tmp2 = pd.Series(list(y))
                mean1 = round(np.mean(list(x)), 3)
                mean2 = round(np.mean(list(y)), 3)
                median1 = round(np.quantile(list(x), 0.5), 3)
                median2 = round(np.quantile(list(y), 0.5), 3)
                fig, (ax1, ax2) = plt.subplots(2,1, figsize=(8, 10))
                tmp = pd.DataFrame({'with {}\nN={}\nmean: {}\nmedian: {}'.format(mut_n, str(len(x)), str(mean1), str(median1)): tmp1,
                                    'without {}\nN={}\nmean: {}\nmedian: {}'.format('mut', str(len(y)), str(mean2), str(median2)): tmp2})

                ax1.set_ylabel('log2(mut/wt)', fontsize=15)
                ax1.set_title('mRNA fold change (with 0)',fontsize=12)
                ax1.tick_params(axis='x', labelsize=14)
                sns.boxplot(data=tmp, ax=ax1, showfliers=False, width=0.5, linewidth=1, showmeans=0,medianprops=dict(color='orange'),
                            order=tmp.columns.to_list()[::-1], whis=1.5, boxprops = {'facecolor':'pink', 'alpha':0.7}, whiskerprops={'linestyle':'--'})
                for p in ax1.artists:
                    b, o, g, a = p.get_facecolor()
                    p.set_facecolor((b, o, g, 0.3))


                out = U_test(list(x), list(y))
                U_m = np.format_float_scientific(out[1], precision = 1)
                U_c = np.format_float_scientific(out[2], precision = 1)
                out = T_test(list(x), list(y))
                T_m = np.format_float_scientific(out[1], precision = 1)
                T_c = np.format_float_scientific(out[2], precision = 1)
                out = KS_test(list(x), list(y))
                KS_m = np.format_float_scientific(out[1], precision = 1)
                KS_c = np.format_float_scientific(out[2], precision = 1)

                text = 'U test: {} > {}: {}'.format(mut_n, 'com',U_m)+'\nU test: {} < {}: {}'.format(mut_n, 'com', U_c,)+'\n-----------------------------------------\nT test: {} > {}: {}'.format(mut_n, 'com', T_m)+'\nT test: {} < {}: {}'.format(mut_n, 'com', T_c)+'\n-----------------------------------------\nKS test: {} > {}: {}'.format(mut_n, 'com', KS_m)+'\nKS test: {} < {}: {}'.format(mut_n, 'com', KS_c)
                #print(text)
                ax1.text(1.6,0,text,fontsize=14, verticalalignment='baseline')


                x = [n for n in d1[abu+'_without0'] if n != 'NULL']
                y = [n for n in d2[abu+'_without0'] if n != 'NULL']
                tmp1 = pd.Series(list(x))
                tmp2 = pd.Series(list(y))
                if len(x)==0:
                    mean1 = 'N/A'
                    median1 = 'N/A'
                else:
                    mean1 = round(np.mean(list(x)), 3)
                    median1 = round(np.quantile(list(x), 0.5), 3)
                if len(y)==0:
                    mean2 = 'N/A'
                    median2 = 'N/A'
                else:
                    mean2 = round(np.mean(list(y)), 3)
                    median2 = round(np.quantile(list(y), 0.5), 3)
                tmp = pd.DataFrame({'with {}\nN={}\nmean: {}\nmedian: {}'.format(mut_n, str(len(x)), str(mean1), str(median1)): tmp1,
                                    'without {}\nN={}\nmean: {}\nmedian: {}'.format('mut', str(len(y)), str(mean2), str(median2)): tmp2})

                ax2.set_ylabel('log2(mut/wt)', fontsize=15)
                ax2.set_title('mRNA fold change (without 0)',fontsize=12)
                ax2.tick_params(axis='x', labelsize=14)
                sns.boxplot(data=tmp, ax=ax2, showfliers=False, width=0.5, linewidth=1, showmeans=0,medianprops=dict(color='orange'),
                            order=tmp.columns.to_list()[::-1], whis=1.5, boxprops = {'facecolor':'pink', 'alpha':0.7}, whiskerprops={'linestyle':'--'})
                for p in ax2.artists:
                    b, o, g, a = p.get_facecolor()
                    p.set_facecolor((b, o, g, 0.3))


                out = U_test(list(x), list(y))
                U_m = np.format_float_scientific(out[1], precision = 1)
                U_c = np.format_float_scientific(out[2], precision = 1)
                out = T_test(list(x), list(y))
                T_m = np.format_float_scientific(out[1], precision = 1)
                T_c = np.format_float_scientific(out[2], precision = 1)
                out = KS_test(list(x), list(y))
                KS_m = np.format_float_scientific(out[1], precision = 1)
                KS_c = np.format_float_scientific(out[2], precision = 1)

                text = 'U test: {} > {}: {}'.format(mut_n, 'com',U_m)+'\nU test: {} < {}: {}'.format(mut_n, 'com', U_c,)+'\n-----------------------------------------\nT test: {} > {}: {}'.format(mut_n, 'com', T_m)+'\nT test: {} < {}: {}'.format(mut_n, 'com', T_c)+'\n-----------------------------------------\nKS test: {} > {}: {}'.format(mut_n, 'com', KS_m)+'\nKS test: {} < {}: {}'.format(mut_n, 'com', KS_c)
                #print(text)
                ax2.text(1.6,0,text,fontsize=14, verticalalignment='baseline')
                plt.tight_layout()
                plt.savefig('figure/abu_plot/FOLD/{}_group{}_{}_{}_{}_with_mutation.{}'.format(d_name, str(group), mut, dn_list[d_n], score_type, fig_type))
                plt.clf()
                plt.close()
                gc.collect()

                d_n += 1

            del d1_list, d2_list
            print()

# mRNA abundance fold change (cdf)
dn_list = ['all data', 'high score', 'middle score', 'low score']

for group in [0]:
    for mut in ['D', 'M']:
        abu_num = 0
        for abu in ['fold_change_avg']:
            abu_num += 1
        
            if mut == 'D':
                mut_n = 'del'
            elif mut == 'M':
                mut_n = 'mis'
            elif mut == 'A':
                mut_n = 'mut'

            ## compare with all mutation (contorol group)
            print('======= compare with all mutation =======')
            # split group
            ana_data = add_two_mRNA_list(data, target, group)
            no_del_clash_result = ana_data[ana_data['A'].astype(str).isin(['[]'])]
            del_clash_result = ana_data[~ana_data[mut].astype(str).isin(['[]'])] 
            ## 無突變資料中將突變資料的 transcript去除 (看有沒有需要)
            no_del_clash_result = no_del_clash_result[~no_del_clash_result['transcript_name'].isin(list(del_clash_result['transcript_name']))]
            ####
            print(len(del_clash_result), len(no_del_clash_result))
            # three different score
            no_one_third_data = no_del_clash_result[no_del_clash_result[score_type]<=one_third]
            no_one_third_data = no_one_third_data[no_one_third_data[score_type]>bot]
            no_two_third_data = no_del_clash_result[no_del_clash_result[score_type]>one_third]
            no_two_third_data = no_two_third_data[no_two_third_data[score_type]<=two_third]
            no_three_third_data = no_del_clash_result[no_del_clash_result[score_type]>two_third]
            no_three_third_data = no_three_third_data[no_three_third_data[score_type]<=top]

            del_one_third_data = del_clash_result[del_clash_result[score_type]>bot]
            del_one_third_data = del_one_third_data[del_one_third_data[score_type]<=one_third]
            del_two_third_data = del_clash_result[del_clash_result[score_type]>one_third]
            del_two_third_data = del_two_third_data[del_two_third_data[score_type]<=two_third]
            del_three_third_data = del_clash_result[del_clash_result[score_type]>two_third]
            del_three_third_data = del_three_third_data[del_three_third_data[score_type]<=top]

            d1_list = [del_clash_result, del_three_third_data, del_two_third_data, del_one_third_data]
            d2_list = [no_del_clash_result, no_three_third_data, no_two_third_data, no_one_third_data]
            d_n = 0
#             print('no {}: '.format('mut'), len(no_three_third_data), len(no_two_third_data), len(no_one_third_data))
#             print('{}: '.format(mut), len(del_three_third_data), len(del_two_third_data), len(del_one_third_data))

            for d1, d2 in zip(d1_list, d2_list):

                print(d_n, title_map_gene[str(group)], dn_list[d_n], mut_n, abu, 'with mutation')
                # plot
                x = d1[abu]
                y = d2[abu]
                tmp1 = pd.Series(list(x))
                tmp2 = pd.Series(list(y))
                
                fig, (ax1, ax2) = plt.subplots(2,1, figsize=(6, 10))
                sorted_x = np.sort(x)
                sorted_y = np.sort(y)
                p_x = 1.*np.arange(len(sorted_x))/float(len(sorted_x)-1)
                p_y = 1.*np.arange(len(sorted_y))/float(len(sorted_y)-1)

                ax1.tick_params(axis='x', labelsize=14)
                ax1.plot(sorted_x, p_x, c='tab:blue', label=mut_n)
                ax1.plot(sorted_y, p_y, c='tab:orange', label='no mutation')
                ax1.legend(loc='upper left', bbox_to_anchor=(1.0, 1.0))
                ax1.set_ylabel('log2(mut/wt)', fontsize=15)
                ax1.set_ylabel('cumulative fraction', fontsize=12)
                ax1.axvline(x=0, c='k', linestyle='dashed', linewidth=0.5)
                out = KS_test(list(x), list(y))
                KS_m = np.format_float_scientific(out[1], precision = 1)
                KS_c = np.format_float_scientific(out[2], precision = 1)

                text = 'KS test: {} > {}: {}'.format(mut_n, 'com', KS_m)+'\nKS test: {} < {}: {}'.format(mut_n, 'com', KS_c)
                #print(text)
                ax1.set_title('mRNA fold change (with 0)\n'+text,fontsize=12)           
                ax1.set_xlim(-1, 1)

                x = [n for n in d1[abu+'_without0'] if n != 'NULL']
                y = [n for n in d2[abu+'_without0'] if n != 'NULL']
                sorted_x = np.sort(x)
                sorted_y = np.sort(y)
                p_x = 1.*np.arange(len(sorted_x))/float(len(sorted_x)-1)
                p_y = 1.*np.arange(len(sorted_y))/float(len(sorted_y)-1)
                ax2.plot(sorted_x, p_x, c='tab:blue', label=mut_n)
                ax2.plot(sorted_y, p_y, c='tab:orange', label='no mutation')
                ax2.legend(loc='upper left', bbox_to_anchor=(1.0, 1.0))
                ax2.set_ylabel('log2(mut/wt)', fontsize=15)
                ax2.tick_params(axis='x', labelsize=14)
                ax2.set_ylabel('cumulative fraction', fontsize=12)
                out = KS_test(list(x), list(y))
                KS_m = np.format_float_scientific(out[1], precision = 1)
                KS_c = np.format_float_scientific(out[2], precision = 1)

                text = 'KS test: {} > {}: {}'.format(mut_n, 'com', KS_m)+'\nKS test: {} < {}: {}'.format(mut_n, 'com', KS_c)
                #print(text)
                ax2.set_title('mRNA fold change (without 0)\n'+text,fontsize=12)
                ax2.axvline(x=0, c='k', linestyle='dashed', linewidth=0.5)
                ax2.set_xlim(-1, 1)
                plt.tight_layout()
                plt.savefig('figure/abu_plot/FOLD_CDF/{}_group{}_{}_{}_{}_with_mutation_cdf.{}'.format(d_name, str(group), mut, d_n, score_type, fig_type))
                plt.clf()
                plt.close()
                gc.collect()

                d_n += 1

            del d1_list, d2_list
            print()

# per score (只分成10等分，一張圖呈現)
score_range = [top, bot]
per_score = (top-bot)/10
#print(per_score)
group_list = [0]
t = data[['transcript_name', 'regulator_name', score_type, 'D', 'M', 'A', 'Gene ID', 'fold_change_avg', 'fold_change_avg_without0']]
top = int(top-per_score)
for group in group_list:
    for mut in ['D', 'M']:
        w = 0
        for w0 in ['fold_change_avg', 'fold_change_avg_without0']:
            if w == 0:
                w_n = 'with0'
                w += 1
            elif w == 1:
                w_n = 'without0'
            print(group, mut, w_n)
            rc_list_mut = []
            rc_list_no = []
            tmp2 = pd.DataFrame()
            my_pal = {}
            p_value = []
            text = []
            text1 = []
            text2 = []
            text3 = []
            order = []
            for s1 in np.arange(top, bot, -int(per_score)):
                s2 = int(s1-per_score)
                tmp = t[t[score_type] <= s1]
                tmp = tmp[tmp[score_type] > s2]
                ana_data = add_two_mRNA_list(tmp, target, group)
                no_del_clash_result = ana_data[ana_data['A'].astype(str).isin(['[]'])]
                del_clash_result = ana_data[~ana_data[mut].astype(str).isin(['[]'])]
                ## 無突變資料中將突變資料的 transcript去除 (看有沒有需要)
                no_del_clash_result = no_del_clash_result[~no_del_clash_result['transcript_name'].isin(list(del_clash_result['transcript_name']))]
                ####
                #x = no_del_clash_result[w0]
                x = pd.Series([n for n in no_del_clash_result[w0] if n != 'NULL'])
                #y = del_clash_result[w0]
                y = pd.Series([n for n in del_clash_result[w0] if n != 'NULL'])
                # test
                out = U_test(list(x), list(y))
                U_m = np.format_float_scientific(out[1], precision = 1)
                U_c = np.format_float_scientific(out[2], precision = 1)
                out = T_test(list(x), list(y))
                T_m = np.format_float_scientific(out[1], precision = 1)
                T_c = np.format_float_scientific(out[2], precision = 1)
                out = KS_test(list(x), list(y))
                KS_m = np.format_float_scientific(out[1], precision = 1)
                KS_c = np.format_float_scientific(out[2], precision = 1)
                #text = '\nU:{}<{}:{}'.format('no', mut, U_c,)+'\n------------------\nT:{}<{}:{}'.format('no', mut, T_c)+'\n------------------\nK:{}<{}:{}'.format('no', mut, KS_c)
                text1.append('U:{}<{}:{}'.format('no', mut, U_c))
                text2.append('T:{}<{}:{}'.format('no', mut, T_c))
                text3.append('K:{}<{}:{}'.format('no', mut, KS_c))
                text.append('U:{}<{}:{}\nU:{}>{}:{}\n---------------------\nT:{}<{}:{}\nT:{}>{}:{}\n---------------------\nK:{}<{}:{}\nK:{}>{}:{}'.format('no', mut, U_c, 'no', mut, U_m,
                                                                                                                'no', mut, T_c, 'no', mut, T_m,                                                                                                              'no', mut, KS_c, 'no', mut, KS_m))
                text1.append('U:{}>{}:{}'.format('no', mut, U_m))
                text2.append('T:{}>{}:{}'.format('no', mut, T_m))
                text3.append('K:{}>{}:{}'.format('no', mut, KS_m))
                
                order.append('No\n                   {}≥s>{}\nN={}'.format(s1, s2, len(x)))
                order.append('Mut\n\nN={}'.format(len(y)))

                tmp3 = pd.DataFrame({'No\n                   {}≥s>{}\nN={}'.format(s1, s2, len(x)): x,
                                     'Mut\n\nN={}'.format(len(y)): y})
                my_pal.update({'No\n                   {}≥s>{}\nN={}'.format(s1, s2, len(x)): 'g',
                                'Mut\n\nN={}'.format(len(y)): 'orange'})
                tmp2 = pd.concat([tmp2, tmp3], axis=1)          

            plt.figure(figsize = (16,8))
            plt.ylabel('mRNA fold change ({})'.format(w_n), fontsize=15)
            plt.tick_params(axis='x', labelsize=11)
            ax = sns.boxplot(data=tmp2, showfliers=False, width=0.5, linewidth=1, order=order,
                             whis=1.5, palette=my_pal, whiskerprops={'linestyle':'--'})
            for p in ax.artists:
                b, o, g, a = p.get_facecolor()
                p.set_facecolor((b, o, g, 0.3))
            columns_list = tmp2.columns

            for i in range(int(len(columns_list)/2)):
                #plt.text(2*(i+1)-1.5,y,text[i],fontsize=11, horizontalalignment='center') 
                plt.annotate(text[i], xy=(i*2/len(columns_list)+0.01,-0.35), xycoords='axes fraction')
            plt.tight_layout()   
            plt.savefig('figure/abu_plot/per_score/{}_group{}_{}_{}_{}_with_mutation.{}'.format(d_name, str(group), mut, score_type, w_n, fig_type))
            plt.clf()
            plt.close()
            gc.collect()
