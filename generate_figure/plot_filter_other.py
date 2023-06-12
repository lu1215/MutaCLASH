from plot_import import *

# per score (只分成10等分，一張圖呈現)

#### 自定義
nor_f = 811.03

# 22G
score_type =  'targeting_score'
top = 10
two_third = 0
one_third = -15
bot = -30

score_range = [top, bot]
per_score = (top-bot)/10
#print(per_score)
group_list = [1,2,8]
t = data[['transcript_name', 'regulator_name', 'targeting_score', '22G_rc_WT', 'D', 'M', 'A', 'Gene ID', 'RNAup_score', '22G_rc_MUT']]
top = int(top-per_score)
for group in group_list:
    for mut in ['D', 'M']:
        print(group, mut)
        rc_list_mut = []
        rc_list_no = []
        tmp2 = pd.DataFrame()
        my_pal = {}
        p_value = []
        text = []
        text1 = []
        text2 = []
        text3 = []
        for s1 in range(top, bot, -int(per_score)):
            s2 = int(s1-per_score)
            #print(group, mut, 'score range: {}-{}'.format(s1, s2))
            tmp = t[t[score_type] <= s1]
            tmp = tmp[tmp[score_type] > s2]
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

            tmp3 = pd.DataFrame({'No\n                  {}≥s>{}\nN={}'.format(s1, s2, len(x)): x,
                                 'Mut\n\nN={}'.format(len(y)): y})
            my_pal.update({'No\n                  {}≥s>{}\nN={}'.format(s1, s2, len(x)): 'g',
                            'Mut\n\nN={}'.format(len(y)): 'orange'})
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
        plt.savefig('figure/G22_plot/per_score/split10/{}_group{}_{}_{}_with_mutation.png'.format(d_name, str(group), mut, score_type))
        plt.clf()
        plt.close()
        gc.collect()

# mRNA
'''
score_type =  'mir_score'
top = max(data[score_type])
two_third = 140
one_third = 100
bot = 60

score_range = [top, bot]
per_score = (top-bot)/10
#print(per_score)
group_list = [0]
t = data[['transcript_name', 'regulator_name', 'mir_score', 'D', 'M', 'A', 'Gene ID', 'fold_change_avg', 'fold_change_avg_without0']]
top = int(top-per_score)
for group in group_list:
    for mut in ['D', 'M']:
        w = 0
        for w0 in ['fold_change_avg', 'fold_change_avg_without0']:
            if w == 0:
                w_n = '_with0'
                w += 1
            elif w == 1:
                w_n = '_without0'
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
            for s1 in range(top, bot, -int(per_score)):
                s2 = int(s1-per_score)
                tmp = t[t[score_type] <= s1]
                tmp = tmp[tmp[score_type] > s2]
                ana_data = add_two_mRNA_list(tmp, group)
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
            plt.savefig('figure/abu_plot/per_score/split10/{}_group{}_{}_{}_{}_with_mutation.png'.format(d_name, str(group), mut, score_type, w_n))
            plt.clf()
            plt.close()
            gc.collect()
'''
