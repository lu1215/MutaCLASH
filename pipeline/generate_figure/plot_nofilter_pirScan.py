from plot_import import *

# all score
tmp_df = data.groupby(by=['targeting_score'])['nor_readcount','nor_count'].sum().reset_index()
for rc_type in ['nor_readcount', 'nor_count']:
    fig, ax = plt.subplots(1,1, figsize=(6, 4), dpi=200)
    ax.plot(tmp_df['targeting_score'], tmp_df[rc_type], color='firebrick')
    ax.set_xlabel('Targeting score')
    ax.set_ylabel('Read counts')
    plt.savefig('figure/{}_{}_pirScan_score.{}'.format(d_name,rc_type, fig_type), bbox_inches='tight')
    plt.clf()
    plt.close()
    gc.collect()

# three part
one_third = np.quantile(list(data['targeting_score']), 0.333)
two_third = np.quantile(list(data['targeting_score']), 0.666)
print("one_third:{}, two_third:{}".format(str(two_third), str(one_third)))
one_third_data = data.copy()
two_third_data = data[data['targeting_score']>=one_third]
three_third_data = data[data['targeting_score']>=two_third]
print("all len:{}, 2/3 len:{}, 1/3 len:{}".format(str(len(one_third_data)), str(len(two_third_data)), str(len(three_third_data))))
print("< Deletion > all len:{}, 2/3 len:{}, 1/3 len:{}".
      format(str(len(one_third_data[~one_third_data['D'].astype(str).isin(['[]'])])), 
      str(len(two_third_data[~two_third_data['D'].astype(str).isin(['[]'])])), 
      str(len(three_third_data[~three_third_data['D'].astype(str).isin(['[]'])]))))
print("< Mismatch > all len:{}, 2/3 len:{}, 1/3 len:{}".
      format(str(len(one_third_data[~one_third_data['M'].astype(str).isin(['[]'])])),
      str(len(two_third_data[~two_third_data['M'].astype(str).isin(['[]'])])),
      str(len(three_third_data[~three_third_data['M'].astype(str).isin(['[]'])]))))

# three part (三分群分開)
'''
one_third = np.quantile(list(data['targeting_score']), 0.333)
two_third = np.quantile(list(data['targeting_score']), 0.666)
print(one_third, two_third)
one_third_data = data[data['targeting_score']<one_third]
two_third_data = data[data['targeting_score']>=one_third]
two_third_data = two_third_data[two_third_data['targeting_score']<two_third]
three_third_data = data[data['targeting_score']>=two_third]
print(len(one_third_data), len(two_third_data), len(three_third_data))
print(len(one_third_data[~one_third_data['D'].astype(str).isin(['[]'])]),
      len(two_third_data[~two_third_data['D'].astype(str).isin(['[]'])]),
      len(three_third_data[~three_third_data['D'].astype(str).isin(['[]'])]))
print(len(one_third_data[~one_third_data['M'].astype(str).isin(['[]'])]),
      len(two_third_data[~two_third_data['M'].astype(str).isin(['[]'])]),
      len(three_third_data[~three_third_data['M'].astype(str).isin(['[]'])]))
'''

# mutation distribution
#### avg ####
# three_third_data、two_third_data、one_third_data
for mut in ['D', 'M']:
    sd = 0
    for split_data in [one_third_data, two_third_data, three_third_data]:
        sd += 1
        for rc_type in ['nor_readcount', 'nor_count']:
            #print(mut, sd, rc_type)
            name_list = []
            rc_list = []
            pi_region = [[] for i in range(61)]
            read_region = [[] for i in range(61)]
            pair_region = [[] for i in range(61)]
            region = [i for i in range(41, -20, -1)]

            input_data = split_data.copy()
            input_data = input_data[~input_data[mut].astype(str).isin(['[]'])]
            input_data.reset_index(drop=True, inplace=True)
            
            rc_mut = 0
            for i in range(len(input_data)):
                rc = input_data[rc_type][i]
                init_pos = int(input_data['pirscan_target_pos'][i].split('-')[1]) + 20
                region_list = input_data['rem_tran_target_pos'][i].split('-')
                start = int(region_list[0])
                stop = int(region_list[1])
                mrns_stop = int(input_data['mRNA_len'][i])
                mut_list = [int(n) for n in eval(input_data[mut][i])]
                for dd in mut_list:
                    rc_mut += rc
                for j in range(61):
                    if mrns_stop <= mrns_stop and init_pos - j >= 0: # 在 mRNA內
                        if init_pos - j > stop:
                            continue
                        elif init_pos - j in mut_list:
                            pi_region[j].append(rc)
                            read_region[j].append(rc)
                            pair_region[j].append(input_data['hybrid_read'][i])
                        elif init_pos - j < start:
                            continue
                        else:
                            read_region[j].append(rc)
                    else:
                        continue

            pi_list = []
            pair_list = []
            total_pi_list = []
            for i in range(61):
                try:
                    total_pi_list.append(sum(pi_region[i]))
                    ave = sum(pi_region[i])/sum(read_region[i]) 
                    pair = len(list(set(pair_region[i]))) 
                except:
                    #print(pi_region[i])
                    ave = 0
                    pair = 0
                pi_list.append(ave)
                pair_list.append(pair)
            total_pi_list = total_pi_list[::-1]
            pi_list = pi_list[::-1]
            pair_list = pair_list[::-1]

            fig, (ax1, ax2) = plt.subplots(2,1, sharex=True, figsize=(6, 10))

            ax1.axis([42, -21, 0, max(pi_list)*1.1])
            ax2.axis([42, -21, 0, max(total_pi_list)*1.1])
            ax1.bar(region, pi_list, color='#FFA07A', width=0.2)
            ax1.scatter(region, pi_list, color='#FFA07A', s=10)
            ax1.set_title('Mutation Distribution'+'\nN='+str(len(input_data)))
            ax1.set_ylabel('Ratio of mutation')
            ax1.axvline(x=1, c='k', linestyle='dashed', linewidth=0.5)
            ax1.axvline(x=21, c='k', linestyle='dashed', linewidth=0.5)
#             ax2.bar(region, pair_list, color='#FFA07A', width=0.2)
#             ax2.scatter(region, pair_list, color='#FFA07A', s=10)
#             ax2.set_title('Mutation Pair Distribution'+'\nN='+str(len(list(set(list(input_data['hybrid_read']))))))
#             ax2.set_ylabel('number of pair')
#             ax2.set_xlabel('pos')
#             ax2.axvline(x=1, c='k', linestyle='dashed', linewidth=0.5)
#             ax2.axvline(x=21, c='k', linestyle='dashed', linewidth=0.5)
            ax2.bar(region, total_pi_list, color='#FFA07A', width=0.2)
            ax2.scatter(region, total_pi_list, color='#FFA07A', s=10)
            ax2.set_title('Mutation Distribution'+'\npercentage in site='+str(round(sum(total_pi_list[20: 41])/rc_mut, 3)*100)+' %')
            ax2.set_ylabel('number of mutation (rc)')
            ax2.set_xlabel('pos')
            ax2.axvline(x=1, c='k', linestyle='dashed', linewidth=0.5)
            ax2.axvline(x=21, c='k', linestyle='dashed', linewidth=0.5)
            if sd == 1:
                sd_type = 'all_data'
            elif sd == 2:  
                sd_type = 'two_third'
            elif sd == 3:
                sd_type = 'one_third'
            plt.savefig('figure/distribution_plot/{}_{}_{}_{}_pirScan.{}'.format(d_name,mut,sd_type,rc_type, fig_type))
            plt.clf()
            plt.close()
            gc.collect()

# pairing ratio

# pos 下有無 deletion 的 acc  (quick)
def pairing_plot(plot_diff=False):
    for mut in ['D', 'M']:
        sd = 0
        for split_data in [one_third_data, two_third_data, three_third_data]:
            sd += 1
            for rc_type in ['nor_readcount', 'nor_count']:
                #print(mut, sd, rc_type)
                d_acc_dict = {}
                acc_dict = {}
                data_acc = split_data.copy()
                reg_seq_len = int(max([len(n) for n in data_acc['regulator_seq']]))
                data_acc.reset_index(drop=True, inplace=True)
                del_data = data_acc[~data_acc[mut].astype(str).isin(['[]'])]
                per_data = data_acc[data_acc[mut].astype(str).isin(['[]'])]
                if plot_diff:
                    ## 20240917 changing original per data(hybrid read with no mutation) to hybrid read with mutation(no matter D or M) ##
                    per_data = data_acc[~data_acc['D'].astype(str).isin(['[]']) | ~data_acc['M'].astype(str).isin(['[]'])]
                    ## ---------------------------------------------------------------------------------------------------------------- ##
                del_data.reset_index(drop=True, inplace=True)
                per_data.reset_index(drop=True, inplace=True)
                #print(len(per_data), len(del_data))

                a_list = []
                t_list = []
                c_list = []
                g_list = []

                match = [0]*reg_seq_len   # 此位置無 deletion的結果
                err = [0]*reg_seq_len
                dash_in_regulator = [0]*reg_seq_len
                dash_in_transcript = [0]*reg_seq_len
                num_in_this_pos = 0
                for i in range(len(per_data)):
                    #print(i, end='\r')
                    num_in_this_pos += 1
                    rc = per_data[rc_type][i] 

                    trans_seq = per_data['pir_target_mRNA_region'][i]
                    reg_seq = per_data['regulator_seq'][i]
                    mir_reg_length = len(trans_seq)
                    shift_dash = 0
                    for p in range(1, mir_reg_length+1):
                        try:
                            if trans_seq[(mir_reg_length-p)].upper() == 'A' and reg_seq[p-1-shift_dash].upper() == 'T':
                                match[p-1-shift_dash] += 1*rc
                                #print(p-1-shift_dash)
                            elif trans_seq[(mir_reg_length-p)].upper() == 'T' and reg_seq[p-1-shift_dash].upper() == 'A': 
                                match[p-1-shift_dash] += 1*rc
                                #print(p-1-shift_dash)
                            elif trans_seq[(mir_reg_length-p)].upper() == 'C' and reg_seq[p-1-shift_dash].upper() == 'G': 
                                match[p-1-shift_dash] += 1*rc
                                #print(p-1-shift_dash)
                            elif trans_seq[(mir_reg_length-p)].upper() == 'G' and reg_seq[p-1-shift_dash].upper() == 'C': 
                                match[p-1-shift_dash] += 1*rc
                                #print(p-1-shift_dash)
                            # GU mismatch
                            elif trans_seq[(mir_reg_length-p)].upper() == 'G' and reg_seq[p-1-shift_dash].upper() == 'T': 
                                match[p-1-shift_dash] += 0.5*rc
                                err[p-1-shift_dash] += 0.5*rc
                            elif trans_seq[(mir_reg_length-p)].upper() == 'T' and reg_seq[p-1-shift_dash].upper() == 'G': 
                                match[p-1-shift_dash] += 0.5*rc
                                err[p-1-shift_dash] += 0.5*rc
                            # else
                            elif trans_seq[(mir_reg_length-p)] == '-':
                                dash_in_transcript[p-1-shift_dash] += 1*rc
                            elif reg_seq[p-1-shift_dash] == '-':
                                dash_in_regulator[p-1-shift_dash] += 1*rc
                                shift_dash += 1
                                reg_seq = reg_seq[:p-shift_dash] + reg_seq[p-shift_dash+1:]

                            else:
                                err[p-1-shift_dash] += 1*rc
                        except:continue
                for pos in range(1, reg_seq_len+1):
                    d_match = [0]*reg_seq_len # 此位置有 deletion的結果
                    d_err = [0]*reg_seq_len

                    dash_in_regulator = [0]*reg_seq_len
                    dash_in_transcript = [0]*reg_seq_len

                    d_num_in_this_pos = 0

                    length = len(del_data)
                    #print(pos, length, end='\r')
                    for i in range(length):

                        init_pos = int(del_data['pirscan_target_pos'][i].split('-')[0])
                        rc = del_data[rc_type][i] 
                        trans_seq = del_data['pir_target_mRNA_region'][i]
                        reg_seq = del_data['regulator_seq'][i]
                        mir_reg_length = len(trans_seq)
                        if del_data[mut][i] != '[]':
                            if len(del_data['regulator_seq'][i]) >= pos: # reg_seq長度要大於等於紀錄的位置
                                D = eval(str(del_data[mut][i]))
                                for d in D:
                                    if (int(d) - init_pos) == (len(del_data['regulator_seq'][i])-pos): # 兩座標反轉
                                        d_num_in_this_pos += 1

                                        for p in range(1, mir_reg_length+1):
                                            try:
                                                #### acc    
                                                if trans_seq[(mir_reg_length-p)].upper() == 'A' and reg_seq[p-1-shift_dash].upper() == 'T':
                                                    d_match[p-1-shift_dash] += 1*rc
                                                    #print(p-1-shift_dash)
                                                elif trans_seq[(mir_reg_length-p)].upper() == 'T' and reg_seq[p-1-shift_dash].upper() == 'A': 
                                                    d_match[p-1-shift_dash] += 1*rc
                                                    #print(p-1-shift_dash)
                                                elif trans_seq[(mir_reg_length-p)].upper() == 'C' and reg_seq[p-1-shift_dash].upper() == 'G': 
                                                    d_match[p-1-shift_dash] += 1*rc
                                                    #print(p-1-shift_dash)
                                                elif trans_seq[(mir_reg_length-p)].upper() == 'G' and reg_seq[p-1-shift_dash].upper() == 'C': 
                                                    d_match[p-1-shift_dash] += 1*rc
                                                    #print(p-1-shift_dash)
                                                # GU mismatch
                                                elif trans_seq[(mir_reg_length-p)].upper() == 'G' and reg_seq[p-1-shift_dash].upper() == 'T': 
                                                    d_match[p-1-shift_dash] += 0.5*rc
                                                    d_err[p-1-shift_dash] += 0.5*rc
                                                elif trans_seq[(mir_reg_length-p)].upper() == 'T' and reg_seq[p-1-shift_dash].upper() == 'G': 
                                                    d_match[p-1-shift_dash] += 0.5*rc
                                                    d_err[p-1-shift_dash] += 0.5*rc

                                                # else
                                                elif trans_seq[(mir_reg_length-p)] == '-':
                                                    dash_in_transcript[p-1-shift_dash] += 1*rc
                                                elif reg_seq[p-1-shift_dash] == '-':
                                                    dash_in_regulator[p-1-shift_dash] += 1*rc
                                                    shift_dash += 1
                                                    reg_seq = reg_seq[:p-shift_dash] + reg_seq[p-shift_dash+1:]
                                                else:
                                                    d_err[p-1-shift_dash] += 1*rc

                                            except:continue
                    d_acc = [d_match[i]/(d_match[i]+d_err[i]) if (d_match[i]+d_err[i])!= 0 else 0 for i in range(reg_seq_len) ]
                    acc = [match[i]/(match[i]+err[i]) if (match[i]+err[i]) != 0 else 0  for i in range(reg_seq_len)] 
                    d_acc_dict.update({str(pos): [d_acc, d_num_in_this_pos]})
                    acc_dict.update({str(pos): [acc, num_in_this_pos]})
                    
                # single plot
                region = [i for i in range(1, reg_seq_len+1)]
                for i in range(1,reg_seq_len+1): 
                    diff = [d_acc_dict[str(i)][0][j]-acc_dict[str(i)][0][j] for j in range(reg_seq_len)]
                    #print(diff[i-1])
                    fig, (ax1, ax2) = plt.subplots(2,1, figsize=(6,8))
                    ax1.axis([reg_seq_len, 1, 0, 1])      
                    ax1.plot(region,acc_dict[str(i)][0], linewidth=0.7, label='pos {}, N={}'.format(str(i), acc_dict[str(i)][1]))
                    ax1.plot(region,d_acc_dict[str(i)][0], linewidth=0.7, label='Mut pos {}, N={}'.format(str(i), d_acc_dict[str(i)][1]), c='orange')
                    ax1.set_xticks(region)  
                    ax1.set_title('mutation (or not) match ratio')
                    ax1.axvline(x=2,c='k',linestyle='dashed', linewidth=0.7)
                    ax1.axvline(x=7,c='k',linestyle='dashed', linewidth=0.7)
                    ax1.axvline(x=i,c='r',linestyle='dashed', linewidth=0.5)
                    ax1.legend(loc='upper left', bbox_to_anchor=(1.0, 1.0))

                    ax2.plot(region, diff, linewidth=0.7, label='diff (mut-all)')
                    ax2.set_title('mutation (mut - all)')
                    ax2.axis([reg_seq_len, 1, min(diff)*1.1, max(diff)*1.1]) 
                    ax2.axvline(x=i,c='r',linestyle='dashed', linewidth=0.5)
                    ax2.axhline(y=0,c='k', linewidth=0.5)
                    ax2.set_xticks(region) 
                    if sd == 1:
                        sd_type = 'all_data'
                    elif sd == 2:  
                        sd_type = 'two_third'
                    elif sd == 3:
                        sd_type = 'one_third'
                    if plot_diff==False:
                        plt.savefig('figure/pairing_ratio_plot/plot21/{}_{}_{}_{}_{}_pirScan.{}'.format(str(i),d_name,mut,sd_type,rc_type, fig_type), bbox_inches='tight')
                    plt.clf()
                    plt.close()
                    gc.collect()
                    
                # pairing ratio V2
                region = [i for i in range(1, reg_seq_len+1)]
                diff_1nt_before = []
                diff_1nt_after = []
                diff_2nt_before = []
                diff_2nt_after = []
                diff_3nt_before = []
                diff_3nt_after = []
                diff_4nt_before = []
                diff_4nt_after = []
                diff_5nt_before = []
                diff_5nt_after = []
                diff_6nt_before = []
                diff_6nt_after = []
                diff_7nt_before = []
                diff_7nt_after = []

                diff_this_pos = []

                for i in range(1,reg_seq_len+1):
                    #print(i, end='\r')
                    if d_acc_dict[str(i)][1] == 0:
                        diff_this_pos.append(0)
                        if i > 1:
                            diff_1nt_before.append(0)
                        if i > 2:
                            diff_2nt_before.append(0)
                        if i > 3:
                            diff_3nt_before.append(0)
                        if i > 4:
                            diff_4nt_before.append(0)
                        if i > 5:
                            diff_5nt_before.append(0)
                        if i > 6:
                            diff_6nt_before.append(0)
                        if i > 7:
                            diff_7nt_before.append(0)

                        if i < reg_seq_len:
                            diff_1nt_after.append(0)
                        if i < reg_seq_len-1:
                            diff_2nt_after.append(0)
                        if i < reg_seq_len-2:
                            diff_3nt_after.append(0)
                        if i < reg_seq_len-3:
                            diff_4nt_after.append(0)
                        if i < reg_seq_len-4:
                            diff_5nt_after.append(0)
                        if i < reg_seq_len-5:
                            diff_6nt_after.append(0)
                        if i < reg_seq_len-6:
                            diff_7nt_after.append(0)
                    else:
                        diff_this_pos.append(d_acc_dict[str(i)][0][i-1]-acc_dict[str(i)][0][i-1])
                        if i > 1:
                            diff_1nt_before.append(d_acc_dict[str(i)][0][i-2]-acc_dict[str(i)][0][i-2])
                        if i > 2:
                            diff_2nt_before.append(d_acc_dict[str(i)][0][i-3]-acc_dict[str(i)][0][i-3])
                        if i > 3:
                            diff_3nt_before.append(d_acc_dict[str(i)][0][i-4]-acc_dict[str(i)][0][i-4])
                        if i > 4:
                            diff_4nt_before.append(d_acc_dict[str(i)][0][i-5]-acc_dict[str(i)][0][i-5])
                        if i > 5:
                            diff_5nt_before.append(d_acc_dict[str(i)][0][i-6]-acc_dict[str(i)][0][i-6])
                        if i > 6:
                            diff_6nt_before.append(d_acc_dict[str(i)][0][i-7]-acc_dict[str(i)][0][i-7])
                        if i > 7:
                            diff_7nt_before.append(d_acc_dict[str(i)][0][i-8]-acc_dict[str(i)][0][i-8])

                        if i < reg_seq_len:
                            diff_1nt_after.append(d_acc_dict[str(i)][0][i]-acc_dict[str(i)][0][i])
                        if i < reg_seq_len-1:
                            diff_2nt_after.append(d_acc_dict[str(i)][0][i+1]-acc_dict[str(i)][0][i+1])
                        if i < reg_seq_len-2:
                            diff_3nt_after.append(d_acc_dict[str(i)][0][i+2]-acc_dict[str(i)][0][i+2])
                        if i < reg_seq_len-3:
                            diff_4nt_after.append(d_acc_dict[str(i)][0][i+3]-acc_dict[str(i)][0][i+3])
                        if i < reg_seq_len-4:
                            diff_5nt_after.append(d_acc_dict[str(i)][0][i+4]-acc_dict[str(i)][0][i+4])
                        if i < reg_seq_len-5:
                            diff_6nt_after.append(d_acc_dict[str(i)][0][i+5]-acc_dict[str(i)][0][i+5])
                        if i < reg_seq_len-6:
                            diff_7nt_after.append(d_acc_dict[str(i)][0][i+6]-acc_dict[str(i)][0][i+6])
                            
                # box
                my_pal = {'minus 7':'g', 'minus 6':'g','minus 5':'g', 'minus 4':'g',
                        'minus 3':'g', 'minus 2':'g','minus 1':'g', 'current\npos':'b',
                        'plus 1':'g', 'plus 2':'g', 'plus 3':'g','plus 4':'g', 'plus 5':'g','plus 6':'g', 'plus 7':'g'}

                tmp2 = pd.DataFrame({'plus 7':pd.Series(diff_7nt_after),'plus 6':pd.Series(diff_6nt_after),
                                    'plus 5':pd.Series(diff_5nt_after),'plus 4':pd.Series(diff_4nt_after),
                                    'plus 3':pd.Series(diff_3nt_after),'plus 2':pd.Series(diff_2nt_after),
                                    'plus 1':pd.Series(diff_1nt_after),'current\npos':pd.Series(diff_this_pos),
                                    'minus 1':pd.Series(diff_1nt_before),'minus 2':pd.Series(diff_2nt_before),
                                    'minus 3':pd.Series(diff_3nt_before), 'minus 4':pd.Series(diff_4nt_before),
                                    'minus 5':pd.Series(diff_5nt_before),'minus 6':pd.Series(diff_6nt_before),
                                    'minus 7':pd.Series(diff_7nt_before)})
                plt.figure(figsize = (15,6))
                ax = sns.boxplot(data=tmp2,
                                order=['plus 7','plus 6','plus 5','plus 4','plus 3','plus 2','plus 1','current\npos',
                                        'minus 1','minus 2','minus 3','minus 4','minus 5','minus 6','minus 7',],
                                showfliers = False, width=0.5, linewidth=1, showmeans=True, palette=my_pal) 
                ax.set_ylabel('Difference', fontsize=14)
                ax.axhline(y=0,c='k', linewidth=0.7, linestyle='--')
                for p in ax.artists:
                    b, o, g, a = p.get_facecolor()
                    p.set_facecolor((b, o, g, 0.3))
                    
                if sd == 1:
                    sd_type = 'all_data'
                elif sd == 2:  
                    sd_type = 'two_third'
                elif sd == 3:
                    sd_type = 'one_third'
                if plot_diff==False:
                    plt.savefig('figure/pairing_ratio_plot/plot_all/{}_{}_{}_{}_pirScan.{}'.format(d_name,mut,sd_type,rc_type, fig_type))
                plt.clf()
                plt.close()
                gc.collect()

                # pairing ratio V3
                for box in ['seed', 'non-seed']:
                    for diff in [True, False]:
                        pos = 21 # maximum x-axis position
                        new_df = pd.DataFrame()
                        for i in range(pos+1):
                            if diff:
                                if i==0:
                                    continue
                                else:
                                    y = [d_acc_dict[str(i)][0][j]-acc_dict[str(i)][0][j] for j in range(reg_seq_len)]
                            else:
                                if i==0:
                                    y = acc_dict['1'][0]
                                else:
                                    y = d_acc_dict[str(i)][0]
                            tmp_df = pd.DataFrame()
                            tmp_df['x'] = [i]*reg_seq_len
                            tmp_df['y'] = y
                            tmp_df['region'] = ['non-seed']+['seed']*6+['non-seed']*(reg_seq_len-7)
                            new_df = pd.concat([new_df, tmp_df])
                        plot_df = new_df[new_df['region'] == box]

                        fig, ax = plt.subplots(1,1, figsize=(15,4), dpi=200)
                        if diff:
                            palette = ['#9abbff']*pos
                        else:
                            palette = ['#dddcdc'] + ['#9abbff']*pos
                            
                        ax = sns.boxplot(data=plot_df, x='x', y='y', width=0.3, showfliers=False,
                                        showmeans=False, medianprops=dict(color='orange'), palette=palette)
                        ax.set_ylabel('pairing ratio')
                        ax.set_xlabel('mutation position')

                        if sd == 1:
                            sd_type = 'all_data'
                        elif sd == 2:  
                            sd_type = 'two_third'
                        elif sd == 3:
                            sd_type = 'one_third'

                        if diff:
                            ## 20240918 add horizontal line to show median of leftmost bar ##
                            med = np.median(plot_df[plot_df['x'] == 0]['y'].to_numpy())
                            ax.axhline(med, linestyle='--', label='avg')
                            ## ----------------------------------------------------------- ##
                            ax.set_title(box)
                            if plot_diff:
                                plt.savefig('figure/pairing_ratio_plot/plot_diff/{}_{}_{}_{}_{}_diff_pirScan.{}'.format(d_name, mut, sd_type, rc_type, box, fig_type), bbox_inches='tight')
                        else:
                            ## 20240918 add horizontal line to show median of leftmost bar ##
                            med = np.median(plot_df[plot_df['x'] == 0]['y'].to_numpy())
                            ax.axhline(med, linestyle='--', label='avg')
                            ## ----------------------------------------------------------- ##
                            ax.set_title(box, y=1.07)
                            ax.axvline(x=0.5,c='k',linestyle='dashed', linewidth=0.7)
                            labels = [item.get_text() for item in ax.get_xticklabels()]
                            labels[0] = 'ALL'
                            ax.set_xticklabels(labels)

                            # statistical test (U-test)
                            ylim = ax.get_ylim()
                            test_pos = ylim[1] + 0.02 * (ylim[1] - ylim[0])
                            for i in range(1, pos+1):
                                sample1 = plot_df.loc[ plot_df['x']==0, 'y']
                                sample2 = plot_df.loc[ plot_df['x']==i, 'y']
                                _, p = stats.mannwhitneyu(sample1, sample2)
                                
                                if p < 0.05:
                                    ax.text(i-0.15, test_pos, '**')
                                elif p < 0.1:
                                    ax.text(i-0.05, test_pos, '*')
                            if plot_diff:
                                plt.savefig('figure/pairing_ratio_plot/plot_diff/{}_{}_{}_{}_{}_pirScan.{}'.format(d_name, mut, sd_type, rc_type, box, fig_type), bbox_inches='tight')
                        
                        plt.clf()
                        plt.close()
                        gc.collect()

pairing_plot(plot_diff=False)
pairing_plot(plot_diff=True)