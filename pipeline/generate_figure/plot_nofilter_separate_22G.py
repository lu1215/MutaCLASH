from plot_import import *

# target_col, count, targeting_score, mir_score
target_col = 'read_count'

score_type = target_col
top = max(data[score_type])
bot = 0

# three part
one_third = np.quantile(list(data[target_col]), 0.333)
two_third = np.quantile(list(data[target_col]), 0.666)
print("one_third:{}, two_third:{}".format(str(two_third), str(one_third)))
one_third_data = data.copy()
two_third_data = data[data[target_col]>=one_third]
three_third_data = data[data[target_col]>=two_third]
print("all len:{}, 2/3 len:{}, 1/3 len:{}".format(str(len(one_third_data)), str(len(two_third_data)), str(len(three_third_data))))
print("< Deletion > all len:{}, 2/3 len:{}, 1/3 len:{}".
      format(str(len(one_third_data[~one_third_data['D'].astype(str).isin(['[]'])])), 
      str(len(two_third_data[~two_third_data['D'].astype(str).isin(['[]'])])), 
      str(len(three_third_data[~three_third_data['D'].astype(str).isin(['[]'])]))))
print("< Mismatch > all len:{}, 2/3 len:{}, 1/3 len:{}".
      format(str(len(one_third_data[~one_third_data['M'].astype(str).isin(['[]'])])),
      str(len(two_third_data[~two_third_data['M'].astype(str).isin(['[]'])])),
      str(len(three_third_data[~three_third_data['M'].astype(str).isin(['[]'])]))))


print('high: {} < S <= {}'.format(str(two_third), str(top)))
print('mid: {} < S <= {}'.format(str(one_third), str(two_third)))
print('low: {} < S <= {}'.format(str(bot), str(one_third)))

# fold change 欄位
alpha1 = find_alpha(data['22G_rc_WT'])
alpha2 = find_alpha(data['22G_rc_MUT'])
alpha = min([alpha1, alpha2])
print(alpha)
data['22G_rc_WT_alpha'] = [data['22G_rc_WT'][i]+alpha for i in range(len(data))]
data['22G_rc_MUT_alpha'] = [data['22G_rc_MUT'][i]+alpha for i in range(len(data))]
data['fold_change'] = [math.log(data['22G_rc_MUT_alpha'][i]/data['22G_rc_WT_alpha'][i], 2) for i in range(len(data))]
print(len(data))
x1 = [n for n in data['22G_rc_MUT']]
y1 = [n for n in data['22G_rc_WT']]
no_0_list = []
for i,j in zip(x1, y1):
    if i != 0 and j != 0:
        no_0_list.append(math.log(i/j, 2))
    else:
        no_0_list.append("NULL")
print(len(no_0_list))
data['fold_change_without0'] = no_0_list

# 22G read count 分布

dn_list = ['all data', 'high score', 'middle score', 'low score']

for group in [1,2,8]:
    for mut in ['D','M']:
        if mut == 'D':
            mut_n = 'del'
        elif mut == 'M':
            mut_n = 'mis'
        elif mut == 'A':
            mut_n = 'mut'
            
#         # compare with deletion/mismatch (control group)
#         print('======= compare with all deletion/mismatch =======')
#         # split group
#         ana_data = add_two_mRNA_list(data, group)
#         no_del_clash_result = ana_data[ana_data[mut].astype(str).isin(['[]'])]
#         del_clash_result = ana_data[~ana_data[mut].astype(str).isin(['[]'])] 
#         print(len(del_clash_result), len(no_del_clash_result))
#         # three different score
#         no_one_third_data = no_del_clash_result[no_del_clash_result[score_type]<=one_third]
#         no_one_third_data = no_one_third_data[no_one_third_data[score_type]>bot]
#         no_two_third_data = no_del_clash_result[no_del_clash_result[score_type]>one_third]
#         no_two_third_data = no_two_third_data[no_two_third_data[score_type]<=two_third]
#         no_three_third_data = no_del_clash_result[no_del_clash_result[score_type]>two_third]
#         no_three_third_data = no_three_third_data[no_three_third_data[score_type]<=top]

#         del_one_third_data = del_clash_result[del_clash_result[score_type]>bot]
#         del_one_third_data = del_one_third_data[del_one_third_data[score_type]<=one_third]
#         del_two_third_data = del_clash_result[del_clash_result[score_type]>one_third]
#         del_two_third_data = del_two_third_data[del_two_third_data[score_type]<=two_third]
#         del_three_third_data = del_clash_result[del_clash_result[score_type]>two_third]
#         del_three_third_data = del_three_third_data[del_three_third_data[score_type]<=top]
        
#         d1_list = [del_clash_result, del_three_third_data, del_two_third_data, del_one_third_data]
#         d2_list = [no_del_clash_result, no_three_third_data, no_two_third_data, no_one_third_data]
#         d_n = 0
#         print('no {}: '.format(mut), len(no_three_third_data), len(no_two_third_data), len(no_one_third_data))
#         print('{}: '.format(mut), len(del_three_third_data), len(del_two_third_data), len(del_one_third_data))
        
#         for d1, d2 in zip(d1_list, d2_list):
    
#             print(d_n, title_map_gene[str(group)], dn_list[d_n], mut_n)
#             # plot
#             x = [math.log((n/nor_f)*1000000, 10) if n != 0 else 0 for n in d1['22G_rc_WT']]
#             y = [math.log((n/nor_f)*1000000, 10) if n != 0 else 0 for n in d2['22G_rc_WT']]
#             tmp1 = pd.Series(list(x))
#             tmp2 = pd.Series(list(y))
#             mean1 = round(np.mean(list(x)), 3)
#             mean2 = round(np.mean(list(y)), 3)
#             median1 = round(np.quantile(list(x), 0.5), 3)
#             median2 = round(np.quantile(list(y), 0.5), 3)
#             tmp = pd.DataFrame({'with {}\nN={}\nmean: {}\nmedian: {}'.format(mut_n, str(len(x)), str(mean1), str(median1)): tmp1,
#                                 'without {}\nN={}\nmean: {}\nmedian: {}'.format(mut_n, str(len(y)), str(mean2), str(median2)): tmp2})

#             plt.figure(figsize = (6,5))
#             plt.ylabel('22G read count', fontsize=15)
#             plt.tick_params(axis='x', labelsize=14)
#             ax = sns.boxplot(data=tmp, showfliers=False, width=0.5, linewidth=1, showmeans=0,medianprops=dict(color='orange'),
#                              whis=1.5, boxprops = {'color':'black', 'facecolor':'pink', 'alpha':0.7}, whiskerprops={'linestyle':'--'})
#             # add_stat_annotation(ax,data=tmp,
#             #                    box_pairs=[('with\nmut','without\nmut')],
#             #                    test='Wilcoxon test', text_format='full', loc='outside', verbose=2,fontsize=12)
#             for p in ax.artists:
#                 b, o, g, a = p.get_facecolor()
#                 p.set_facecolor((b, o, g, 0.3))


#             out = U_test(list(x), list(y))
#             U_m = np.format_float_scientific(out[1], precision = 1)
#             U_c = np.format_float_scientific(out[2], precision = 1)
#             out = T_test(list(x), list(y))
#             T_m = np.format_float_scientific(out[1], precision = 1)
#             T_c = np.format_float_scientific(out[2], precision = 1)
#             out = KS_test(list(x), list(y))
#             KS_m = np.format_float_scientific(out[1], precision = 1)
#             KS_c = np.format_float_scientific(out[2], precision = 1)

#             text = 'U test: {} > {}: {}'.format(mut_n, 'com',U_m)+'\nU test: {} < {}: {}'.format(mut_n, 'com', U_c,)+'\n-----------------------------------------\nT test: {} > {}: {}'.format(mut_n, 'com', T_m)+'\nT test: {} < {}: {}'.format(mut_n, 'com', T_c)+'\n-----------------------------------------\nKS test: {} > {}: {}'.format(mut_n, 'com', KS_m)+'\nKS test: {} < {}: {}'.format(mut_n, 'com', KS_c)
#             #print(text)
#             plt.text(1.6,0,text,fontsize=14)
#             plt.savefig('figure/G22_plot/22G/{}_group{}_{}_{}_{}.{}'.format(d_name, str(group), mut, d_n, score_type, fig_type), bbox_inches='tight')
#             plt.clf()
#             plt.close()
#             gc.collect()
            
#             d_n += 1
            
#         del d1_list, d2_list
#         print()

        
        ## compare with all mutation (contorol group)
        print('======= compare with all mutation =======')
        # split group
        ana_data = add_two_mRNA_list(data, group)
        no_del_clash_result = ana_data[ana_data['A'].astype(str).isin(['[]'])]
        del_clash_result = ana_data[~ana_data[mut].astype(str).isin(['[]'])] 
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
        print('no {}: '.format('mut'), len(no_three_third_data), len(no_two_third_data), len(no_one_third_data))
        print('{}: '.format(mut), len(del_three_third_data), len(del_two_third_data), len(del_one_third_data))
        
        for d1, d2 in zip(d1_list, d2_list):
    
            print(d_n, title_map_gene[str(group)], dn_list[d_n], mut_n, 'with mutation')
            # plot
#             x = [math.log((n/nor_f)*1000000, 10) if n != 0 else 0 for n in d1['22G_rc_WT']]
#             y = [math.log((n/nor_f)*1000000, 10) if n != 0 else 0 for n in d2['22G_rc_WT']]
            x = [n/nor_f for n in d1['22G_rc_WT']]
            y = [n/nor_f for n in d2['22G_rc_WT']]
            tmp1 = pd.Series(list(x))
            tmp2 = pd.Series(list(y))
            mean1 = round(np.mean(list(x)), 3)
            mean2 = round(np.mean(list(y)), 3)
            median1 = round(np.quantile(list(x), 0.5), 3)
            median2 = round(np.quantile(list(y), 0.5), 3)
            tmp = pd.DataFrame({'with {}\nN={}\nmean: {}\nmedian: {}'.format(mut_n, str(len(x)), str(mean1), str(median1)): tmp1,
                                'without {}\nN={}\nmean: {}\nmedian: {}'.format('mut', str(len(y)), str(mean2), str(median2)): tmp2})

            plt.figure(figsize = (6,5))
            plt.ylabel('22G read count', fontsize=15)
            plt.tick_params(axis='x', labelsize=14)
            ax = sns.boxplot(data=tmp, showfliers=False, width=0.5, linewidth=1, showmeans=0,medianprops=dict(color='orange'),
                             order=tmp.columns.to_list()[::-1], whis=1.5, boxprops = {'facecolor':'pink', 'alpha':0.7}, whiskerprops={'linestyle':'--'})
            # add_stat_annotation(ax,data=tmp,
            #                    box_pairs=[('with\nmut','without\nmut')],
            #                    test='Wilcoxon test', text_format='full', loc='outside', verbose=2,fontsize=12)
            for p in ax.artists:
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
            plt.text(1.6,0,text,fontsize=14)
            plt.savefig('figure/G22_plot/22G/{}_group{}_{}_{}_{}_with_mutation.{}'.format(d_name, str(group), mut, d_n, score_type, fig_type), bbox_inches='tight')
            plt.clf()
            plt.close()
            gc.collect()
            
            d_n += 1
            
        del d1_list, d2_list
        print()

# FOLD CHANGE

dn_list = ['all data', 'high score', 'middle score', 'low score']

for group in [1,2,8]:
    for mut in ['D', 'M']:
        if mut == 'D':
            mut_n = 'del'
        elif mut == 'M':
            mut_n = 'mis'
        elif mut == 'A':
            mut_n = 'mut'
            
#         # compare with deletion/mismatch (control group)
#         print('======= compare with all deletion/mismatch =======')
#         # split group
#         ana_data = add_two_mRNA_list(data, group)
#         no_del_clash_result = ana_data[ana_data[mut].astype(str).isin(['[]'])]
#         del_clash_result = ana_data[~ana_data[mut].astype(str).isin(['[]'])] 
#         print(len(del_clash_result), len(no_del_clash_result))
#         # three different score
#         no_one_third_data = no_del_clash_result[no_del_clash_result[score_type]<=one_third]
#         no_one_third_data = no_one_third_data[no_one_third_data[score_type]>bot]
#         no_two_third_data = no_del_clash_result[no_del_clash_result[score_type]>one_third]
#         no_two_third_data = no_two_third_data[no_two_third_data[score_type]<=two_third]
#         no_three_third_data = no_del_clash_result[no_del_clash_result[score_type]>two_third]
#         no_three_third_data = no_three_third_data[no_three_third_data[score_type]<=top]

#         del_one_third_data = del_clash_result[del_clash_result[score_type]>bot]
#         del_one_third_data = del_one_third_data[del_one_third_data[score_type]<=one_third]
#         del_two_third_data = del_clash_result[del_clash_result[score_type]>one_third]
#         del_two_third_data = del_two_third_data[del_two_third_data[score_type]<=two_third]
#         del_three_third_data = del_clash_result[del_clash_result[score_type]>two_third]
#         del_three_third_data = del_three_third_data[del_three_third_data[score_type]<=top]
        
#         d1_list = [del_clash_result, del_three_third_data, del_two_third_data, del_one_third_data]
#         d2_list = [no_del_clash_result, no_three_third_data, no_two_third_data, no_one_third_data]
#         d_n = 0
#         print('no {}: '.format(mut), len(no_three_third_data), len(no_two_third_data), len(no_one_third_data))
#         print('{}: '.format(mut), len(del_three_third_data), len(del_two_third_data), len(del_one_third_data))
        
#         for d1, d2 in zip(d1_list, d2_list):
    
#             print(d_n, title_map_gene[str(group)], dn_list[d_n], mut_n)
#             # plot
#             x = d1['fold_change']
#             y = d2['fold_change']
#             tmp1 = pd.Series(list(x))
#             tmp2 = pd.Series(list(y))
#             mean1 = round(np.mean(list(x)), 3)
#             mean2 = round(np.mean(list(y)), 3)
#             median1 = round(np.quantile(list(x), 0.5), 3)
#             median2 = round(np.quantile(list(y), 0.5), 3)
#             tmp = pd.DataFrame({'with {}\nN={}\nmean: {}\nmedian: {}'.format(mut_n, str(len(x)), str(mean1), str(median1)): tmp1,
#                                 'without {}\nN={}\nmean: {}\nmedian: {}'.format(mut_n, str(len(y)), str(mean2), str(median2)): tmp2})
            
#             fig, (ax1, ax2) = plt.subplots(2,1, figsize=(8, 10))
#             ax1.set_ylabel('log2 (wt+α/mut+α)', fontsize=15)
#             ax1.tick_params(axis='x', labelsize=14)
#             sns.boxplot(data=tmp, ax=ax1, showfliers=False, width=0.5, linewidth=1, showmeans=0,medianprops=dict(color='orange'),
#                              whis=1.5, boxprops = {'color':'black', 'facecolor':'pink', 'alpha':0.7}, whiskerprops={'linestyle':'--'})
#             for p in ax1.artists:
#                 b, o, g, a = p.get_facecolor()
#                 p.set_facecolor((b, o, g, 0.3))


#             out = U_test(list(x), list(y))
#             U_m = np.format_float_scientific(out[1], precision = 1)
#             U_c = np.format_float_scientific(out[2], precision = 1)
#             out = T_test(list(x), list(y))
#             T_m = np.format_float_scientific(out[1], precision = 1)
#             T_c = np.format_float_scientific(out[2], precision = 1)
#             out = KS_test(list(x), list(y))
#             KS_m = np.format_float_scientific(out[1], precision = 1)
#             KS_c = np.format_float_scientific(out[2], precision = 1)
            
#             text = 'U test: {} > {}: {}'.format(mut_n, 'com',U_m)+'\nU test: {} < {}: {}'.format(mut_n, 'com', U_c,)+'\n-----------------------------------------\nT test: {} > {}: {}'.format(mut_n, 'com', T_m)+'\nT test: {} < {}: {}'.format(mut_n, 'com', T_c)+'\n-----------------------------------------\nKS test: {} > {}: {}'.format(mut_n, 'com', KS_m)+'\nKS test: {} < {}: {}'.format(mut_n, 'com', KS_c)
#             #print(text)
#             ax1.text(1.6,0,text,fontsize=14, verticalalignment='top')
            
            
#             x = [n for n in d1['fold_change_without0'] if n != 'NULL']
#             y = [n for n in d2['fold_change_without0'] if n != 'NULL']
#             tmp1 = pd.Series(list(x))
#             tmp2 = pd.Series(list(y))
#             mean1 = round(np.mean(list(x)), 3)
#             mean2 = round(np.mean(list(y)), 3)
#             median1 = round(np.quantile(list(x), 0.5), 3)
#             median2 = round(np.quantile(list(y), 0.5), 3)
#             tmp = pd.DataFrame({'with {}\nN={}\nmean: {}\nmedian: {}'.format(mut_n, str(len(x)), str(mean1), str(median1)): tmp1,
#                                 'without {}\nN={}\nmean: {}\nmedian: {}'.format(mut_n, str(len(y)), str(mean2), str(median2)): tmp2})

#             ax2.set_ylabel('log2 (mut/wt)', fontsize=15)
#             ax2.tick_params(axis='x', labelsize=14)
#             sns.boxplot(data=tmp, ax=ax2, showfliers=False, width=0.5, linewidth=1, showmeans=0,medianprops=dict(color='orange'),
#                              whis=1.5, boxprops = {'color':'black', 'facecolor':'pink', 'alpha':0.7}, whiskerprops={'linestyle':'--'})
#             for p in ax2.artists:
#                 b, o, g, a = p.get_facecolor()
#                 p.set_facecolor((b, o, g, 0.3))


#             out = U_test(list(x), list(y))
#             U_m = np.format_float_scientific(out[1], precision = 1)
#             U_c = np.format_float_scientific(out[2], precision = 1)
#             out = T_test(list(x), list(y))
#             T_m = np.format_float_scientific(out[1], precision = 1)
#             T_c = np.format_float_scientific(out[2], precision = 1)
#             out = KS_test(list(x), list(y))
#             KS_m = np.format_float_scientific(out[1], precision = 1)
#             KS_c = np.format_float_scientific(out[2], precision = 1)

#             text = 'U test: {} > {}: {}'.format(mut_n, 'com',U_m)+'\nU test: {} < {}: {}'.format(mut_n, 'com', U_c,)+'\n-----------------------------------------\nT test: {} > {}: {}'.format(mut_n, 'com', T_m)+'\nT test: {} < {}: {}'.format(mut_n, 'com', T_c)+'\n-----------------------------------------\nKS test: {} > {}: {}'.format(mut_n, 'com', KS_m)+'\nKS test: {} < {}: {}'.format(mut_n, 'com', KS_c)
#             #print(text)
#             ax2.text(1.6,0,text,fontsize=14, verticalalignment='top')
#             plt.tight_layout()
#             plt.savefig('figure/G22_plot/22G/group{}_{}_{}_{}_without0.{}'.format(str(group), mut, d_n, score_type, fig_type))
#             plt.clf()
#             plt.close()
#             gc.collect()
            
#             d_n += 1
            
#         del d1_list, d2_list
#         print()

        
        ## compare with all mutation (contorol group)
        print('======= compare with all mutation =======')
        # split group
        ana_data = add_two_mRNA_list(data, group)
        no_del_clash_result = ana_data[ana_data['A'].astype(str).isin(['[]'])]
        del_clash_result = ana_data[~ana_data[mut].astype(str).isin(['[]'])] 
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
        print('no {}: '.format('mut'), len(no_three_third_data), len(no_two_third_data), len(no_one_third_data))
        print('{}: '.format(mut), len(del_three_third_data), len(del_two_third_data), len(del_one_third_data))
        
        for d1, d2 in zip(d1_list, d2_list):
    
            print(d_n, title_map_gene[str(group)], dn_list[d_n], mut_n, 'with mutation')
            # plot
            x = d1['fold_change']
            y = d2['fold_change']
            tmp1 = pd.Series(list(x))
            tmp2 = pd.Series(list(y))
            mean1 = round(np.mean(list(x)), 3)
            mean2 = round(np.mean(list(y)), 3)
            median1 = round(np.quantile(list(x), 0.5), 3)
            median2 = round(np.quantile(list(y), 0.5), 3)
            fig, (ax1, ax2) = plt.subplots(2,1, figsize=(8, 10))
            tmp = pd.DataFrame({'with {}\nN={}\nmean: {}\nmedian: {}'.format(mut_n, str(len(x)), str(mean1), str(median1)): tmp1,
                                'without {}\nN={}\nmean: {}\nmedian: {}'.format('mut', str(len(y)), str(mean2), str(median2)): tmp2})
            
            ax1.set_ylabel('log2 (mut+α/wt+α)', fontsize=15)
            ax1.set_title('22G FOLD Change (with 0)',fontsize=12)
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
            ax1.text(1.6,0,text,fontsize=14, verticalalignment='top')
            
            
            x = [n for n in d1['fold_change_without0'] if n != 'NULL']
            y = [n for n in d2['fold_change_without0'] if n != 'NULL']
            tmp1 = pd.Series(list(x))
            tmp2 = pd.Series(list(y))
            mean1 = round(np.mean(list(x)), 3)
            mean2 = round(np.mean(list(y)), 3)
            median1 = round(np.quantile(list(x), 0.5), 3)
            median2 = round(np.quantile(list(y), 0.5), 3)
            tmp = pd.DataFrame({'with {}\nN={}\nmean: {}\nmedian: {}'.format(mut_n, str(len(x)), str(mean1), str(median1)): tmp1,
                                'without {}\nN={}\nmean: {}\nmedian: {}'.format('mut', str(len(y)), str(mean2), str(median2)): tmp2})
        
            ax2.set_ylabel('log2 (mut/wt)', fontsize=15)
            ax2.set_title('22G FOLD Change (without 0)',fontsize=12)
            ax2.tick_params(axis='x', labelsize=14)
            sns.boxplot(data=tmp, ax=ax2, showfliers=False, width=0.5, linewidth=1, showmeans=0,medianprops=dict(color='orange'),
                             whis=1.5, boxprops = {'facecolor':'pink', 'alpha':0.7}, whiskerprops={'linestyle':'--'})
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
            ax2.text(1.6,0,text,fontsize=14, verticalalignment='top')
            plt.tight_layout()
            plt.savefig('figure/G22_plot/FOLD/{}_group{}_{}_{}_{}_with_mutation.{}'.format(d_name, str(group), mut, d_n, score_type, fig_type))
            plt.clf()
            plt.close()
            gc.collect()
            
            d_n += 1
            
        del d1_list, d2_list
        print()

# 22G CDF

dn_list = ['all data', 'high score', 'middle score', 'low score']
#g22_rc_filter = 10
for group in [1,2,8]:
    for mut in ['D','M']:
        if mut == 'D':
            mut_n = 'del'
        elif mut == 'M':
            mut_n = 'mis'
        elif mut == 'A':
            mut_n = 'mut'
        data_t = data.copy()
#         # compare with deletion/mismatch (control group)
#         print('======= compare with all deletion/mismatch =======')
#         # split group

#         ana_data = add_two_mRNA_list(data_t, group)
#         no_del_clash_result = ana_data[ana_data[mut].astype(str).isin(['[]'])]
#         del_clash_result = ana_data[~ana_data[mut].astype(str).isin(['[]'])] 
#         print(len(del_clash_result), len(no_del_clash_result))
#         # three different score
#         no_one_third_data = no_del_clash_result[no_del_clash_result[score_type]<=one_third]
#         no_one_third_data = no_one_third_data[no_one_third_data[score_type]>bot]
#         no_two_third_data = no_del_clash_result[no_del_clash_result[score_type]>one_third]
#         no_two_third_data = no_two_third_data[no_two_third_data[score_type]<=two_third]
#         no_three_third_data = no_del_clash_result[no_del_clash_result[score_type]>two_third]
#         no_three_third_data = no_three_third_data[no_three_third_data[score_type]<=top]

#         del_one_third_data = del_clash_result[del_clash_result[score_type]>bot]
#         del_one_third_data = del_one_third_data[del_one_third_data[score_type]<=one_third]
#         del_two_third_data = del_clash_result[del_clash_result[score_type]>one_third]
#         del_two_third_data = del_two_third_data[del_two_third_data[score_type]<=two_third]
#         del_three_third_data = del_clash_result[del_clash_result[score_type]>two_third]
#         del_three_third_data = del_three_third_data[del_three_third_data[score_type]<=top]
        
#         d1_list = [del_clash_result, del_three_third_data, del_two_third_data, del_one_third_data]
#         d2_list = [no_del_clash_result, no_three_third_data, no_two_third_data, no_one_third_data]
#         d_n = 0
#         print('no {}: '.format(mut), len(no_three_third_data), len(no_two_third_data), len(no_one_third_data))
#         print('{}: '.format(mut), len(del_three_third_data), len(del_two_third_data), len(del_one_third_data))
        
#         for d1, d2 in zip(d1_list, d2_list):
#             x = [math.log((n/nor_f)*1000000, 10) if n != 0 else 0 for n in d1['22G_rc_WT']]
#             y = [math.log((n/nor_f)*1000000, 10) if n != 0 else 0 for n in d2['22G_rc_WT']]
#             print(len(x), len(y))
#             sorted_x = np.sort(x)
#             sorted_y = np.sort(y)
#             p_x = 1.*np.arange(len(sorted_x))/float(len(sorted_x)-1)
#             p_y = 1.*np.arange(len(sorted_y))/float(len(sorted_y)-1)
#             if mut == 'D':
#                 mut_n = 'deletion'
#                 no_mut_n = 'no deletion'
#             elif mut == 'M':
#                 mut_n = 'mismatch'
#                 no_mut_n = 'no mismatch'

#             fig, (ax1, ax2) = plt.subplots(2,1, figsize=(6, 8))
#             #ax1.axis([-1, g22_rc_filter*1.1, 0, 1])
#             ax1.plot(sorted_x, p_x, c='tab:blue', label=mut_n)
#             ax1.plot(sorted_y, p_y, c='tab:orange', label=no_mut_n)
#             ax1.legend(loc='upper left', bbox_to_anchor=(1.0, 1.0))
#             ax1.set_xlabel('log10 (22G read count)', fontsize=12)
#             ax1.set_ylabel('cumulative fraction', fontsize=12)
#             out = KS_test(list(d1['22G_rc_WT']), list(d2['22G_rc_WT']))
#             KS_m = np.format_float_scientific(out[1], precision = 1)
#             KS_c = np.format_float_scientific(out[2], precision = 1)

#             text = 'KS test: {} > {}: {}'.format(mut_n, 'com', KS_m)+'\nKS test: {} < {}: {}'.format(mut_n, 'com', KS_c)
#             #print(text)
#             ax1.set_title(text,fontsize=12)
            

#             x = d1['targeting_score']
#             y = d2['targeting_score']
#             print(len(x), len(y))
#             sorted_x = np.sort(x)
#             sorted_y = np.sort(y)
#             p_x = 1.*np.arange(len(sorted_x))/float(len(sorted_x)-1)
#             p_y = 1.*np.arange(len(sorted_y))/float(len(sorted_y)-1)

#             ax2.plot(sorted_x, p_x, c='tab:blue', label=mut_n)
#             ax2.plot(sorted_y, p_y, c='tab:orange', label=no_mut_n)
#             ax2.legend(loc='upper left', bbox_to_anchor=(1.0, 1.0))
#             ax2.set_xlabel('targeting score', fontsize=12)
#             ax2.set_ylabel('cumulative fraction', fontsize=12)
#             print(d_n, title_map_gene[str(group)], dn_list[d_n], mut_n)
            
#             # plot
#             out = KS_test(list(x), list(y))
#             KS_m = np.format_float_scientific(out[1], precision = 1)
#             KS_c = np.format_float_scientific(out[2], precision = 1)
            
#             text = 'KS test: {} > {}: {}'.format(mut_n, 'com', KS_m)+'\nKS test: {} < {}: {}'.format(mut_n, 'com', KS_c)
#             #print(text)
#             ax2.set_title(text,fontsize=12)
#             plt.tight_layout()
#             plt.savefig('figure/G22_plot/tmp/{}_group{}_{}_{}_{}.{}'.format(d_name, str(group), mut, d_n, score_type, fig_type))
#             plt.clf()
#             plt.close()
#             gc.collect()
            
#             d_n += 1
            
#         del d1_list, d2_list
#         print()

        
        ## compare with all mutation (contorol group)
        print('======= compare with all mutation =======')
        # split group
        ana_data = add_two_mRNA_list(data_t, group)
        no_del_clash_result = ana_data[ana_data['A'].astype(str).isin(['[]'])]
        del_clash_result = ana_data[~ana_data[mut].astype(str).isin(['[]'])] 
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
        print('no {}: '.format('mut'), len(no_three_third_data), len(no_two_third_data), len(no_one_third_data))
        print('{}: '.format(mut), len(del_three_third_data), len(del_two_third_data), len(del_one_third_data))
        
        for d1, d2 in zip(d1_list, d2_list):
    
            print(d_n, title_map_gene[str(group)], dn_list[d_n], mut_n, 'with mutation')
            # plot
            x = [math.log((n/nor_f)*1000000, 10) if n != 0 else 0 for n in d1['22G_rc_WT']]
            y = [math.log((n/nor_f)*1000000, 10) if n != 0 else 0 for n in d2['22G_rc_WT']]
            print(len(x), len(y))
            sorted_x = np.sort(x)
            sorted_y = np.sort(y)
            p_x = 1.*np.arange(len(sorted_x))/float(len(sorted_x)-1)
            p_y = 1.*np.arange(len(sorted_y))/float(len(sorted_y)-1)
            if mut == 'D':
                mut_n = 'deletion'
                no_mut_n = 'no deletion'
            elif mut == 'M':
                mut_n = 'mismatch'
                no_mut_n = 'no mismatch'

            fig, (ax1, ax2) = plt.subplots(2,1, figsize=(6, 8))
            #ax1.axis([-1, g22_rc_filter*1.1, 0, 1])
            ax1.plot(sorted_x, p_x, c='tab:blue', label=mut_n)
            ax1.plot(sorted_y, p_y, c='tab:orange', label='no mutation')
            ax1.legend(loc='upper left', bbox_to_anchor=(1.0, 1.0))
            
            ax1.set_xlabel('log10 (22G read count)', fontsize=12)
            ax1.set_ylabel('cumulative fraction', fontsize=12)
            out = KS_test(list(d1['22G_rc_WT']), list(d2['22G_rc_WT']))
            KS_m = np.format_float_scientific(out[1], precision = 1)
            KS_c = np.format_float_scientific(out[2], precision = 1)

            text = 'KS test: {} > {}: {}'.format(mut_n, 'com', KS_m)+'\nKS test: {} < {}: {}'.format(mut_n, 'com', KS_c)
            #print(text)
            ax1.set_title('22G read count CDF\n'+text, fontsize=12)
            
            # target score cdf
            x = d1['targeting_score']
            y = d2['targeting_score']
            print(len(x), len(y))
            sorted_x = np.sort(x)
            sorted_y = np.sort(y)
            p_x = 1.*np.arange(len(sorted_x))/float(len(sorted_x)-1)
            p_y = 1.*np.arange(len(sorted_y))/float(len(sorted_y)-1)

            ax2.plot(sorted_x, p_x, c='tab:blue', label=mut_n)
            ax2.plot(sorted_y, p_y, c='tab:orange', label='no mutation')
            ax2.legend(loc='upper left', bbox_to_anchor=(1.0, 1.0))
            ax2.set_xlabel('targeting score', fontsize=12)
            ax2.set_ylabel('cumulative fraction', fontsize=12)
            print(d_n, title_map_gene[str(group)], dn_list[d_n], mut_n)
            out = KS_test(list(x), list(y))
            KS_m = np.format_float_scientific(out[1], precision = 1)
            KS_c = np.format_float_scientific(out[2], precision = 1)
            
            text = 'KS test: {} > {}: {}'.format(mut_n, 'com', KS_m)+'\nKS test: {} < {}: {}'.format(mut_n, 'com', KS_c)
            #print(text)
            ax2.set_title('targeting score CDF\n'+text,fontsize=12)
            plt.tight_layout()
            
            plt.savefig('figure/G22_plot/22G_CDF/{}_group{}_{}_{}_{}_with_mutation.{}'.format(d_name, str(group), mut, d_n, score_type, fig_type))
            plt.clf()
            plt.close()
            gc.collect()
            
            d_n += 1
            
        del d1_list, d2_list
        print()

# FOLD CHANGE (cdf)

dn_list = ['all data', 'high score', 'middle score', 'low score']

for group in [1,2,8]:
    for mut in ['D','M']:
        if mut == 'D':
            mut_n = 'del'
        elif mut == 'M':
            mut_n = 'mis'
        elif mut == 'A':
            mut_n = 'mut'
            
#         # compare with deletion/mismatch (control group)
#         print('======= compare with all deletion/mismatch =======')
#         # split group
#         ana_data = add_two_mRNA_list(data, group)
#         no_del_clash_result = ana_data[ana_data[mut].astype(str).isin(['[]'])]
#         del_clash_result = ana_data[~ana_data[mut].astype(str).isin(['[]'])] 
#         print(len(del_clash_result), len(no_del_clash_result))
#         # three different score
#         no_one_third_data = no_del_clash_result[no_del_clash_result[score_type]<=one_third]
#         no_one_third_data = no_one_third_data[no_one_third_data[score_type]>bot]
#         no_two_third_data = no_del_clash_result[no_del_clash_result[score_type]>one_third]
#         no_two_third_data = no_two_third_data[no_two_third_data[score_type]<=two_third]
#         no_three_third_data = no_del_clash_result[no_del_clash_result[score_type]>two_third]
#         no_three_third_data = no_three_third_data[no_three_third_data[score_type]<=top]

#         del_one_third_data = del_clash_result[del_clash_result[score_type]>bot]
#         del_one_third_data = del_one_third_data[del_one_third_data[score_type]<=one_third]
#         del_two_third_data = del_clash_result[del_clash_result[score_type]>one_third]
#         del_two_third_data = del_two_third_data[del_two_third_data[score_type]<=two_third]
#         del_three_third_data = del_clash_result[del_clash_result[score_type]>two_third]
#         del_three_third_data = del_three_third_data[del_three_third_data[score_type]<=top]
        
#         d1_list = [del_clash_result, del_three_third_data, del_two_third_data, del_one_third_data]
#         d2_list = [no_del_clash_result, no_three_third_data, no_two_third_data, no_one_third_data]
#         d_n = 0
#         print('no {}: '.format(mut), len(no_three_third_data), len(no_two_third_data), len(no_one_third_data))
#         print('{}: '.format(mut), len(del_three_third_data), len(del_two_third_data), len(del_one_third_data))
        
#         for d1, d2 in zip(d1_list, d2_list):
    
#             print(d_n, title_map_gene[str(group)], dn_list[d_n], mut_n)
#             # plot
#             x = d1['fold_change']
#             y = d2['fold_change']
#             sorted_x = np.sort(x)
#             sorted_y = np.sort(y)
#             p_x = 1.*np.arange(len(sorted_x))/float(len(sorted_x)-1)
#             p_y = 1.*np.arange(len(sorted_y))/float(len(sorted_y)-1)
            
#             fig, (ax1, ax2) = plt.subplots(2,1, figsize=(8, 10))
            
#             ax1.tick_params(axis='x', labelsize=14)
            
#             ax1.plot(sorted_x, p_x, c='tab:blue', label=mut_n)
#             ax1.plot(sorted_y, p_y, c='tab:orange', label='no mutation')
#             ax1.legend(loc='upper left', bbox_to_anchor=(1.0, 1.0))
#             ax1.set_ylabel('log2 (wt+α/mut+α)', fontsize=15)
#             ax1.set_ylabel('cumulative fraction', fontsize=12)
#             out = KS_test(list(x), list(y))
#             KS_m = np.format_float_scientific(out[1], precision = 1)
#             KS_c = np.format_float_scientific(out[2], precision = 1)
            
#             text = 'KS test: {} > {}: {}'.format(mut_n, 'com', KS_m)+'\nKS test: {} < {}: {}'.format(mut_n, 'com', KS_c)
#             #print(text)
#             ax1.set_title(text,fontsize=12)   
#             ax1.axvline(x=0, c='k', linestyle='dashed', linewidth=0.5)
            
#             x = [n for n in d1['fold_change_without0'] if n != 'NULL']
#             y = [n for n in d2['fold_change_without0'] if n != 'NULL']
            
#             sorted_x = np.sort(x)
#             sorted_y = np.sort(y)
#             p_x = 1.*np.arange(len(sorted_x))/float(len(sorted_x)-1)
#             p_y = 1.*np.arange(len(sorted_y))/float(len(sorted_y)-1)

#             ax2.plot(sorted_x, p_x, c='tab:blue', label=mut_n)
#             ax2.plot(sorted_y, p_y, c='tab:orange', label='no mutation')
#             ax2.legend(loc='upper left', bbox_to_anchor=(1.0, 1.0))
#             ax2.set_ylabel('log2 (mut/wt)', fontsize=15)
#             ax2.tick_params(axis='x', labelsize=14)
#             ax2.set_ylabel('cumulative fraction', fontsize=12)
#             out = KS_test(list(x), list(y))
#             KS_m = np.format_float_scientific(out[1], precision = 1)
#             KS_c = np.format_float_scientific(out[2], precision = 1)
            
#             text = 'KS test: {} > {}: {}'.format(mut_n, 'com', KS_m)+'\nKS test: {} < {}: {}'.format(mut_n, 'com', KS_c)
#             #print(text)
#             ax2.set_title(text,fontsize=12)
#             ax2.axvline(x=0, c='k', linestyle='dashed', linewidth=0.5)
#             plt.tight_layout()
#             plt.savefig('figure/G22_plot/FOLD_CDF/{}_group{}_{}_{}_{}_without0.{}'.format(d_name, str(group), mut, d_n, score_type, fig_type))
#             plt.clf()
#             plt.close()
#             gc.collect()
            
#             d_n += 1
            
#         del d1_list, d2_list
#         print()

        
        ## compare with all mutation (contorol group)
        print('======= compare with all mutation =======')
        # split group
        ana_data = add_two_mRNA_list(data, group)
        no_del_clash_result = ana_data[ana_data['A'].astype(str).isin(['[]'])]
        del_clash_result = ana_data[~ana_data[mut].astype(str).isin(['[]'])] 
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
        print('no {}: '.format('mut'), len(no_three_third_data), len(no_two_third_data), len(no_one_third_data))
        print('{}: '.format(mut), len(del_three_third_data), len(del_two_third_data), len(del_one_third_data))
        
        for d1, d2 in zip(d1_list, d2_list):
    
            print(d_n, title_map_gene[str(group)], dn_list[d_n], mut_n, 'with mutation')
            # plot
            x = d1['fold_change']
            y = d2['fold_change']
            tmp1 = pd.Series(list(x))
            tmp2 = pd.Series(list(y))
            mean1 = round(np.mean(list(x)), 3)
            mean2 = round(np.mean(list(y)), 3)
            median1 = round(np.quantile(list(x), 0.5), 3)
            median2 = round(np.quantile(list(y), 0.5), 3)
            fig, (ax1, ax2) = plt.subplots(2,1, figsize=(8, 10))
            sorted_x = np.sort(x)
            sorted_y = np.sort(y)
            p_x = 1.*np.arange(len(sorted_x))/float(len(sorted_x)-1)
            p_y = 1.*np.arange(len(sorted_y))/float(len(sorted_y)-1)
            
            fig, (ax1, ax2) = plt.subplots(2,1, figsize=(8, 10))
            
            ax1.tick_params(axis='x', labelsize=14)
            ax1.plot(sorted_x, p_x, c='tab:blue', label=mut_n)
            ax1.plot(sorted_y, p_y, c='tab:orange', label='no mutation')
            ax1.legend(loc='upper left', bbox_to_anchor=(1.0, 1.0))
            ax1.set_ylabel('log2 (mut+α/wt+α)', fontsize=15)
            ax1.set_ylabel('cumulative fraction', fontsize=12)
            ax1.axvline(x=0, c='k', linestyle='dashed', linewidth=0.5)
            out = KS_test(list(x), list(y))
            KS_m = np.format_float_scientific(out[1], precision = 1)
            KS_c = np.format_float_scientific(out[2], precision = 1)
            
            text = 'KS test: {} > {}: {}'.format(mut_n, 'com', KS_m)+'\nKS test: {} < {}: {}'.format(mut_n, 'com', KS_c)
            #print(text)
            ax1.set_title('22G FOLD Change (with 0)\n'+text,fontsize=12)           
            
            
            x = [n for n in d1['fold_change_without0'] if n != 'NULL']
            y = [n for n in d2['fold_change_without0'] if n != 'NULL']
            sorted_x = np.sort(x)
            sorted_y = np.sort(y)
            p_x = 1.*np.arange(len(sorted_x))/float(len(sorted_x)-1)
            p_y = 1.*np.arange(len(sorted_y))/float(len(sorted_y)-1)
            ax2.plot(sorted_x, p_x, c='tab:blue', label=mut_n)
            ax2.plot(sorted_y, p_y, c='tab:orange', label='no mutation')
            ax2.legend(loc='upper left', bbox_to_anchor=(1.0, 1.0))
            ax2.set_ylabel('log2 (mut/wt)', fontsize=15)
            ax2.tick_params(axis='x', labelsize=14)
            ax2.set_ylabel('cumulative fraction', fontsize=12)
            out = KS_test(list(x), list(y))
            KS_m = np.format_float_scientific(out[1], precision = 1)
            KS_c = np.format_float_scientific(out[2], precision = 1)
            
            text = 'KS test: {} > {}: {}'.format(mut_n, 'com', KS_m)+'\nKS test: {} < {}: {}'.format(mut_n, 'com', KS_c)
            #print(text)
            ax2.set_title('22G FOLD Change (without 0)\n'+text,fontsize=12)
            ax2.axvline(x=0, c='k', linestyle='dashed', linewidth=0.5)
            plt.tight_layout()
            plt.savefig('figure/G22_plot/FOLD_CDF/{}_group{}_{}_{}_{}_with_mutation.{}'.format(d_name, str(group), mut, d_n, score_type, fig_type))
            plt.clf()
            plt.close()
            gc.collect()
            
            d_n += 1
            
        del d1_list, d2_list

# 每五分數區間，目前只有 pirScan targeting score
'''
# score_range = [5, 0] # [5, 0], [0, -5], [-5, -10], [-10, -15], [-15, -20], [-20, -25], [-25, -30]
# group = 1
# mut = 'M'
t = data[['transcript_name', 'regulator_name', 'targeting_score', '22G_rc_WT', 'D', 'M', 'A', 'Gene ID', 'RNAup_score', '22G_rc_MUT']]
score_list = [[5, 0], [0, -5], [-5, -10], [-10, -15], [-15, -20], [-20, -25], [-25, -30]]
group_list = [1,2,8]
for group in group_list:
    for score_range in score_list:
        for mut in ['D', 'M']:
            print(group, mut, score_range)
            rc_list_mut = []
            rc_list_no = []
            tmp2 = pd.DataFrame()
            my_pal = {}
            p_value = []
            text = []
            text1 = []
            text2 = []
            text3 = []
            for s1 in range(score_range[0], score_range[1], -1):
                for s2 in range(2):
                    ss = s1 - 0.5*s2
                    tmp = t[t['targeting_score'] == ss]
                    tmp['22G_rc_WT'] = [n/nor_f for n in tmp['22G_rc_WT']]
                    ana_data = add_two_mRNA_list(tmp, group)
                    no_del_clash_result = ana_data[ana_data['A'].astype(str).isin(['[]'])]
                    del_clash_result = ana_data[~ana_data[mut].astype(str).isin(['[]'])]
                    x = no_del_clash_result['22G_rc_WT']
                    y = del_clash_result['22G_rc_WT']
                    
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
                                                                                                                    'no', mut, T_c, 'no', mut, T_m,
                                                                                                                   'no', mut, KS_c, 'no', mut, KS_m))
                    text1.append('U:{}>{}:{}'.format('no', mut, U_m))
                    text2.append('T:{}>{}:{}'.format('no', mut, T_m))
                    text3.append('K:{}>{}:{}'.format('no', mut, KS_m))

                    tmp3 = pd.DataFrame({'No\n{}\nN={}'.format(ss, len(x)): x,
                                         'Mut\n{}\nN={}'.format(ss, len(y)): y})
                    my_pal.update({'No\n{}\nN={}'.format(ss, len(x)): 'g',
                                    'Mut\n{}\nN={}'.format(ss, len(y)): 'orange'})
                    tmp2 = pd.concat([tmp2, tmp3], axis=1)

            plt.figure(figsize = (16,8))
            plt.ylabel('22G read count', fontsize=15)
            plt.tick_params(axis='x', labelsize=11)
            ax = sns.boxplot(data=tmp2, showfliers=False, width=0.5, linewidth=1, 
                             whis=1.5, palette=my_pal, whiskerprops={'linestyle':'--'})
            for p in ax.artists:
                b, o, g, a = p.get_facecolor()
                p.set_facecolor((b, o, g, 0.3))
            columns_list = tmp2.columns

            for i in range(int(len(columns_list)/2)):
                #plt.text(2*(i+1)-1.5,y,text[i],fontsize=11, horizontalalignment='center') 
                plt.annotate(text[i], xy=(i*2/len(columns_list)+0.01,-0.35), xycoords='axes fraction')
            plt.tight_layout()   
            plt.savefig('figure/G22_plot/per_score/{}_group{}_{}_{}_{}_with_mutation.{}'.format(d_name, str(group), mut, score_type, str(score_range), fig_type))
            plt.clf()
            plt.close()
            gc.collect()
'''

