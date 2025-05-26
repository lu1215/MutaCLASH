import pandas as pd

def filter_hybrid_seq(df: pd.DataFrame) -> pd.DataFrame:
    """
    按 hybrid_seq 群組，依據 D 與 M 的空/非空組態篩選列。
    條件見 docstring 上方說明。
    """
    def _filter_group(g: pd.DataFrame) -> pd.DataFrame:
        both_empty_mask = (g["D"] == "[]") & (g["M"] == "[]")
        
        # 三種情況判斷
        if both_empty_mask.all():                       # ① 全部都空 → 保留
            return g
        elif (~both_empty_mask).all():                  # ② 全部至少有一欄不空 → 保留
            return g
        else:                                           # ③ 混合 → 只留空的那幾列
            return g[both_empty_mask]
    
    # group_keys=False 可避免把 group 標籤變成 MultiIndex
    return df.groupby("hybrid_seq", group_keys=False).apply(_filter_group)