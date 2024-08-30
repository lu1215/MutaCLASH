import argparse
import os
from tqdm import tqdm, trange


file_path = os.path.dirname(os.path.abspath(__file__))

def process_dup_seq_and_count(input_path: str, data_name: str):
    with open(input_path, 'r') as f:
        store_next_line = False
        sequence_dict = {}
        for line in tqdm(f, desc="Processing input sequences"):
            if line.startswith('@'):
                store_next_line = True
            elif store_next_line:
                line = line.strip('\n')
                if sequence_dict.get(line) == None:
                    sequence_dict[line] = 1
                ## sequence already exists, value plus 1
                else:
                    sequence_dict[line] += 1
                store_next_line = False
    
    output_path = "{}/output/{}.fa".format(file_path, data_name)
    with open(output_path, 'w') as f_out:
        key_list = sorted(list(sequence_dict.keys()))
        print("Processing output:\t")
        pbar = tqdm(total=len(key_list))
        for idx, seq in zip(range(len(sequence_dict.keys())), key_list):
            f_out.write(">{}_{}\n{}\n".format(idx, sequence_dict[seq], seq))
            pbar.update(1)
    # print("number of unique reads: {}".format(len(sequence_dict.keys())))      
    # print("read counts: {}".format(sum(sequence_dict.values())))      

# def main():

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    # parser.add_argument("--read", help="read file path", type=str)
    # parser.add_argument("--regulator", help="regulator file path", type=str)
    # parser.add_argument("--target", help="target file path", type=str)
    # # parser.add_argument("--tool", help="tool chira/clan/hyb", type=str)
    # parser.add_argument("--tool", help="tool (only chira now)", type=str)
    parser.add_argument("--data_path", help="input data_path", type=str)
    parser.add_argument("--data_name", help="input data_name", type=str)
    args = parser.parse_args()
    inputpath = args.data_path
    inputname = args.data_name
    # outputname = args.output
    # typename = args.type
    # parser = argparse.ArgumentParser()
    # parser.add_argument("--input", help="input file path", type=str)
    # parser.add_argument("--output", help="output file path", type=str)
    # args = parser.parse_args()
    # input_path = args.input
    # input_path = "/media/disk1/shangyi/mutaclash_0723_test/ALG-1_rep3_trimmed.fq"
    # output_path = args.output
    process_dup_seq_and_count(inputpath, inputname)
    
