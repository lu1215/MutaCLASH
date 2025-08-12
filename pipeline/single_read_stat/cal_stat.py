import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu



def construct_rc_table(df_target, df_single_D_rc, df_single_M_rc,
                       tname_col_target='Gene name', tname_col_reads='transcript_name',
                       pos_col='pos', count_col='readcount', pos_is_1_based=True):
    # 初始化兩個 coverage array 欄位
    df_target = df_target.copy()
    df_target['D_rc_arr'] = df_target['length'].apply(lambda l: np.zeros(l, dtype=float))
    df_target['M_rc_arr'] = df_target['length'].apply(lambda l: np.zeros(l, dtype=float))

    # 建立 transcript -> df row index 映射
    name_to_row = pd.Series(df_target.index.to_numpy(), index=df_target[tname_col_target]).to_dict()

    def accumulate(df_reads, dest_col):
        # 把 reads 表轉成 numpy，並做越界過濾
        codes = df_reads[tname_col_reads].map(lambda x: name_to_row.get(x, -1)).to_numpy()
        pos0 = df_reads[pos_col].to_numpy(np.int64) - (1 if pos_is_1_based else 0)
        counts = df_reads[count_col].to_numpy(float)

        lengths = df_target['length'].to_numpy()
        valid = (codes >= 0) & (pos0 >= 0)
        valid &= pos0 < lengths[codes.clip(min=0)]  # 防止 -1 索引

        if not valid.any():
            return

        codes_v  = codes[valid]
        pos_v    = pos0[valid]
        counts_v = counts[valid]

        # 依 df row index 排序 -> 分段 -> 每段一次 np.add.at
        order   = np.argsort(codes_v, kind='stable')
        codes_s = codes_v[order]
        pos_s   = pos_v[order]
        counts_s= counts_v[order]
        cuts = np.flatnonzero(np.r_[True, codes_s[1:] != codes_s[:-1], True])

        for s, e in zip(cuts[:-1], cuts[1:]):
            ridx = codes_s[s]
            arr = df_target.at[ridx, dest_col]
            np.add.at(arr, pos_s[s:e], counts_s[s:e])  # 就地累加

    # 分別累加 D 與 M
    accumulate(df_single_D_rc, 'D_rc_arr')
    accumulate(df_single_M_rc, 'M_rc_arr')
    # print(df_target)

    return df_target

def Utest_pvalue(binding_site, target_name, df_target):
    site_st, site_end = map(int, binding_site.split('-'))
    target_RNA_D_rc = df_target.loc[df_target['Gene name'] == target_name, 'D_rc_arr'].iloc[0]
    target_RNA_M_rc = df_target.loc[df_target['Gene name'] == target_name, 'M_rc_arr'].iloc[0]
    site_D_rc = target_RNA_D_rc[site_st - 1 : site_end]
    site_M_rc = target_RNA_M_rc[site_st - 1 : site_end]
    _, pvalue_D = mannwhitneyu(target_RNA_D_rc, site_D_rc)
    _, pvalue_M = mannwhitneyu(target_RNA_M_rc, site_M_rc)
    return pvalue_D, pvalue_M

def compute_pvalues_for_hybrid(df_target, df_hybrid,
                               site_col='binding_site',
                               baseline='outside'):
    """
    回傳一個新的 df_hybrid（copy）並新增 'pvalue_D', 'pvalue_M' 欄位。
    兼容 Python 3.6：用 .values 取出 numpy 陣列，避免依賴較新 API。
    """
    out = df_hybrid.copy()

    # 建立名稱 -> 陣列的快取（避免每列都用 DataFrame 查找）
    mapD = dict(zip(df_target["Gene name"].values, df_target['D_rc_arr'].values))
    mapM = dict(zip(df_target["Gene name"].values, df_target['M_rc_arr'].values))

    names = out["transcript_name"].values
    sites = out[site_col].values

    pD_list, pM_list = [], []

    def safe_mwu(x, y, alternative='two-sided'):
        x = np.asarray(x)
        y = np.asarray(y)

        # 去除 NaN
        if x.size == 0 or y.size == 0:
            return np.nan
        x = x[~np.isnan(x)]
        y = y[~np.isnan(y)]
        if x.size == 0 or y.size == 0:
            return np.nan

        # 合併後所有值皆相同 -> p=1.0
        allvals = np.concatenate((x, y))
        # np.ptp = max - min
        if np.ptp(allvals) == 0:
            return 1.0

        # 舊版 SciPy 的保險：若仍因 ties 出錯就回 1.0
        try:
            return mannwhitneyu(x, y, alternative=alternative)[1]
        except ValueError as e:
            if 'identical' in str(e):
                return 1.0
            raise
    # mw = mannwhitneyu  # 局部別名，減少全域查找開銷

    for i in range(len(out)):
        tname = names[i]
        bs    = sites[i]

        arrD = mapD.get(tname)
        arrM = mapM.get(tname)

        # 若找不到對應 transcript，記 NaN
        if arrD is None or arrM is None or not isinstance(bs, str) or '-' not in bs:
            pD_list.append(np.nan); pM_list.append(np.nan); continue

        # try:
        st, ed = bs.split('-')
        st = int(st); ed = int(ed)
        # except Exception:
            # pD_list.append(np.nan); pM_list.append(np.nan); continue

        start = max(0, st - 1)
        end   = min(len(arrD), ed)

        siteD = arrD[start:end]
        siteM = arrM[start:end]

        if baseline == 'outside':
            baseD = np.concatenate((arrD[:start], arrD[end:])) if end > start else arrD
            baseM = np.concatenate((arrM[:start], arrM[end:])) if end > start else arrM
        else:
            baseD = arrD
            baseM = arrM

        if siteD.size == 0 or baseD.size == 0:
            pD = np.nan
        else:
            pD = safe_mwu(baseD, siteD, alternative='two-sided')

        if siteM.size == 0 or baseM.size == 0:
            pM = np.nan
        else:
            pM = safe_mwu(baseM, siteM, alternative='two-sided')

        pD_list.append(pD)
        pM_list.append(pM)

    out['pvalue_D'] = np.array(pD_list)
    out['pvalue_M'] = np.array(pM_list)
    return out

def main():
    import argparse
    from os.path import basename, splitext, join
    parser = argparse.ArgumentParser(description="calculating p-value of identified region mutation readcount")
    parser.add_argument("--target", help="target CSV 檔案路徑", required=True)
    parser.add_argument("--hybrid_csv", help="輸入 hybrid metadata CSV 檔案路徑")  # 若未使用可移除
    parser.add_argument("--single_D", help="輸入 single read CSV 檔案路徑", required=True)
    parser.add_argument("--single_M", help="輸入 single read CSV 檔案路徑", required=True)
    parser.add_argument("--output_dir", help="輸出資料夾路徑")
    args = parser.parse_args()

    df_target = pd.read_csv(args.target)
    # 假設有 sequence 欄位
    df_target['length'] = df_target['sequence'].str.len()
    df_hybrid = pd.read_csv(args.hybrid_csv)
    df_single_D_rc = pd.read_csv(args.single_D)
    df_single_M_rc = pd.read_csv(args.single_M)

    # 產生 coverage 陣列（就地回傳）
    df_target = construct_rc_table(df_target, df_single_D_rc, df_single_M_rc)

    # base = splitext(basename(args.single_D))[0]
    # df_target.to_csv(args.output_dir + base.replace("chira_single_step1_detail_D_nor_rc", "with_single_mut_readcount") + ".csv")
    df_out = compute_pvalues_for_hybrid(df_target, df_hybrid, site_col='rem_tran_target_pos', baseline='whole')
    df_out["single_read_mut_stat_significance_D"] = df_out["pvalue_D"] <= 0.01
    df_out["single_read_mut_stat_significance_M"] = df_out["pvalue_M"] <= 0.01
    base = splitext(basename(args.hybrid_csv))[0]
    df_out.to_csv(args.output_dir + base + "_with_sgl_mut_stat.csv", index=False)
    # 之後你就可以：
    # target_RNA_22g = df_target.loc[df_target['Gene name']==target_name, 'WT_HRDE1_22G'].iloc[0]
    # 或者直接用 M/D 陣列做切片
    # site = target_RNA_22g[site_st - 1 : site_end]
    # ...

if __name__ == "__main__":
    main()
