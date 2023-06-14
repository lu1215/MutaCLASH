import argparse
import pandas as pd
import os

parser = argparse.ArgumentParser()
parser.add_argument("--input", help="type input file path", type=str)
parser.add_argument("--reg", help="type regulator file path", type=str)
parser.add_argument("--trans", help="type transcript file path", type=str)
args = parser.parse_args()
input_list = [args.input, args.reg, args.trans]

for i in input_list:
    name, extension = os.path.splitext(i)
    if extension=='.fasta':
        os.system("cp {}.fasta {}.fa".format(name,name))
        os.system
