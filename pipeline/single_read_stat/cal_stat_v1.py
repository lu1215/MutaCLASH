import pandas as pd 
import argparse
import numpy as np
from scipy.stats import mannwhitneyu


def parse_args():
    parser = argparse.ArgumentParser(description="calculating p-value of identified region mutation readcount")
    parser.add_argument("--target", help="target CSV 檔案路徑")
    parser.add_argument("--hybrid_csv", help="輸入 hybrid metadata CSV 檔案路徑")
    parser.add_argument("--single_D", help="輸入 single read  CSV 檔案路徑")
    parser.add_argument("--single_M", help="輸入 hybrid metadata CSV 檔案路徑")
    parser.add_argument("output_dir", help="輸出資料夾路徑")
    return parser.parse_args()

def construct_rc_table(single_read_D_df, single_read_M_df):
    df_target["D_rc_arr"] = df_target["length"].apply(lambda l: np.zeros(l, dtype=float))
    df_target["M_rc_arr"] = df_target["length"].apply(lambda l: np.zeros(l, dtype=float))

    # 向量化相加
    for tname, sub_reads in single_read_D_df.groupby('transcript_name'):
        if tname in df_target['Gene name'].values:
            idx = df_target.index[df_target['Gene name'] == tname][0]
            arr = df_target.at[idx, 'D_rc_array']
            positions = sub_reads['pos'].to_numpy() - 1  # 轉成 0-based
            counts = sub_reads['readcount'].to_numpy()
            np.add.at(arr, positions, counts)  # 就地相加（避免覆蓋）
    
    for tname, sub_reads in single_read_M_df.groupby('transcript_name'):
        if tname in df_target['Gene name'].values:
            idx = df_target.index[df_target['Gene name'] == tname][0]
            arr = df_target.at[idx, 'M_rc_array']
            positions = sub_reads['pos'].to_numpy() - 1  # 轉成 0-based
            counts = sub_reads['readcount'].to_numpy()
            np.add.at(arr, positions, counts)  # 就地相加（避免覆蓋）

    
    

def Utest_pvalue(binding_site, target_name, df_target):
    site_st, site_end = map(int, binding_site.split('-'))
    target_RNA_22g = df_target[df_target['Gene name'] == target_name]['WT_HRDE1_22G'].reset_index(drop=True)[0]
    site_22g = target_RNA_22g[site_st - 1 : site_end]
    statistic, pvalue = mannwhitneyu(target_RNA_22g, site_22g)
    return pvalue

def main():
    args = parse_args()
    global df_target
    df_target = pd.read_csv(args.target)
    df_target["length"] = df_target['sequence'].apply(len)

    df_hybrid_meta = pd.read_csv(args.hybrid_csv)
    df_single_D_rc = pd.read_csv(args.single_D)
    df_single_M_rc = pd.read_csv(args.single_M)




if __name__ == "__main__":
    main()