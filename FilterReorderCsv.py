import pandas as pd
import argparse

def process_csv(input_file, output_file):
    # 讀取 CSV 檔案
    df = pd.read_csv(input_file)
    
    # 保留所需欄位
    selected_columns = [
        "hybrid_seq", "read_count", "D", "M", "transcript_name",
        "regulator_name", "rem_tran_target_pos", "remain_pos", "on_reg_pos", 
        "reg_hyb_target_pos", "targeting_score", "mir_energy", "mir_transcript_seq", 
        "mir_regulator_seq", "RNAup_score", "RNAup_transcript_seq", "RNAup_regulator_seq"
    ]
    df_filtered = df[selected_columns]
    
    # 重新排序
    ordered_columns = [
        "hybrid_seq", "read_count", "D", "M", "transcript_name",
        "regulator_name", "rem_tran_target_pos", "remain_pos", "on_reg_pos", 
        "reg_hyb_target_pos", "targeting_score", "mir_energy", "mir_transcript_seq", 
        "mir_regulator_seq", "RNAup_score", "RNAup_transcript_seq", "RNAup_regulator_seq"
    ]
    df_reordered = df_filtered[ordered_columns]
    
    # 儲存到新的 CSV 檔案
    df_reordered.to_csv(output_file, index=False)
    # print(f"處理完成，輸出檔案：{output_file}")

# 主程式入口
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="處理 CSV 檔案，篩選並重新排序欄位。")
    parser.add_argument("input_csv", type=str, help="輸入的 CSV 檔案名稱")
    parser.add_argument("output_csv", type=str, help="輸出的 CSV 檔案名稱")
    args = parser.parse_args()
    
    process_csv(args.input_csv, args.output_csv)
