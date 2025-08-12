# # -*- coding: utf-8 -*-
# import argparse
# import csv
# import ast
# import os
# from collections import defaultdict

# def parse_args():
#     parser = argparse.ArgumentParser(description="統計 D 與 M 欄位的 nor_readcount")
#     parser.add_argument("input_csv", help="輸入 CSV 檔案路徑")
#     parser.add_argument("output_dir", help="輸出資料夾路徑")
#     return parser.parse_args()

# def safe_parse_list(s):
#     """將字串轉成 list，如果是空的 [] 就回傳空列表"""
#     try:
#         lst = ast.literal_eval(s)
#         if isinstance(lst, list):
#             return lst
#         else:
#             return []
#     except Exception:
#         return []

# def process_file(input_csv, output_dir):
#     d_counts = defaultdict(float)  # (transcript_name, pos) -> nor_readcount 累加
#     m_counts = defaultdict(float)

#     with open(input_csv, "r", newline='', encoding="utf-8") as f:
#         reader = csv.DictReader(f)
#         for row in reader:
#             transcript = row["transcript_name"]
#             nor_readcount = float(row["nor_readcount"])

#             # 處理 D 欄位
#             d_positions = safe_parse_list(row["D"])
#             for pos in d_positions:
#                 if nor_readcount != 0:
#                     d_counts[(transcript, pos)] += nor_readcount

#             # 處理 M 欄位
#             m_positions = safe_parse_list(row["M"])
#             for pos in m_positions:
#                 if nor_readcount != 0:
#                     m_counts[(transcript, pos)] += nor_readcount

#     base_name = os.path.splitext(os.path.basename(input_csv))[0]

#     # 輸出 D 統計
#     d_output_path = os.path.join(output_dir, base_name + "_D_nor_rc.csv")
#     with open(d_output_path, "w", newline='', encoding="utf-8") as f:
#         writer = csv.writer(f)
#         writer.writerow(["transcript_name", "pos", "readcount"])
#         for (transcript, pos), count in sorted(d_counts.items()):
#             writer.writerow([transcript, pos, count])

#     # 輸出 M 統計
#     m_output_path = os.path.join(output_dir, base_name + "_M_nor_rc.csv")
#     with open(m_output_path, "w", newline='', encoding="utf-8") as f:
#         writer = csv.writer(f)
#         writer.writerow(["transcript_name", "pos", "readcount"])
#         for (transcript, pos), count in sorted(m_counts.items()):
#             writer.writerow([transcript, pos, count])

#     print("已輸出：")
#     print(d_output_path)
#     print(m_output_path)

# if __name__ == "__main__":
#     args = parse_args()
#     if not os.path.exists(args.output_dir):
#         os.makedirs(args.output_dir)
#     process_file(args.input_csv, args.output_dir)

# -*- coding: utf-8 -*-
import argparse
import csv
import os
from multiprocessing import Pool, cpu_count
from collections import defaultdict

def parse_args():
    p = argparse.ArgumentParser(description="以多進程統計 D/M 欄位的 nor_readcount")
    p.add_argument("input_csv", help="輸入 CSV 檔案路徑")
    p.add_argument("output_dir", help="輸出資料夾路徑")
    p.add_argument("-w", "--workers", type=int, default=None,
                   help="進程數（未指定時自動使用 CPU 一半，向下取整，至少 1）")
    p.add_argument("-c", "--chunksize", type=int, default=20000,
                   help="每批處理的資料筆數（越大 IPC 開銷越小，但佔用記憶體越多）")
    return p.parse_args()

def fast_parse_positions(s):
    if not s:
        return []
    s = s.strip()
    if s == "[]":
        return []
    if s.startswith("[") and s.endswith("]"):
        inner = s[1:-1].strip()
        if not inner:
            return []
        out = []
        for part in inner.split(","):
            part = part.strip()
            if not part:
                continue
            try:
                out.append(int(part))
            except ValueError:
                try:
                    out.append(int(float(part)))
                except Exception:
                    pass
        return out
    return []

def process_chunk(rows):
    d_part = defaultdict(float)
    m_part = defaultdict(float)
    for transcript, d_s, m_s, nor_s in rows:
        try:
            nor = float(nor_s)
        except Exception:
            continue
        if nor == 0.0:
            continue

        for pos in fast_parse_positions(d_s):
            d_part[(transcript, pos)] += nor
        for pos in fast_parse_positions(m_s):
            m_part[(transcript, pos)] += nor
    return dict(d_part), dict(m_part)

def read_in_chunks(dict_reader, chunk_size):
    chunk = []
    for row in dict_reader:
        chunk.append((row["transcript_name"], row["D"], row["M"], row["nor_readcount"]))
        if len(chunk) >= chunk_size:
            yield chunk
            chunk = []
    if chunk:
        yield chunk

def merge_counts(dst, src):
    for k, v in src.items():
        dst[k] = dst.get(k, 0.0) + v

def write_out(path, counts_dict):
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["transcript_name", "pos", "readcount"])
        for (transcript, pos) in sorted(counts_dict.keys(), key=lambda x: (x[0], x[1])):
            w.writerow([transcript, pos, counts_dict[(transcript, pos)]])

def process_file_mp(input_csv, output_dir, workers, chunksize):
    from os.path import basename, splitext, join
    d_counts, m_counts = {}, {}
    with open(input_csv, "r", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        with Pool(processes=workers) as pool:
            for d_part, m_part in pool.imap_unordered(process_chunk, read_in_chunks(reader, chunksize)):
                merge_counts(d_counts, d_part)
                merge_counts(m_counts, m_part)
    base = splitext(basename(input_csv))[0]
    d_out = join(output_dir, base + "_D_nor_rc.csv")
    m_out = join(output_dir, base + "_M_nor_rc.csv")
    write_out(d_out, d_counts)
    write_out(m_out, m_counts)
    # print("已輸出：\n{}\n{}".format(d_out, m_out))

def main():
    args = parse_args()
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    # 自動偵測：未指定 workers 或 <=0 時，用 CPU 核心數的一半（至少 1）
    if args.workers is None or args.workers <= 0:
        try:
            auto_workers = max(1, int(cpu_count() / 2))
        except Exception:
            auto_workers = 1
        workers = auto_workers
    else:
        workers = args.workers

    process_file_mp(args.input_csv, args.output_dir, workers, args.chunksize)

if __name__ == "__main__":
    main()
