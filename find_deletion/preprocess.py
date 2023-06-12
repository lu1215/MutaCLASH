##################################
# Copyright (C) 2023 Ryan Chung  #
#                                #
# History:                       #
# 2023/01/20 Ryan Chung          #
# 			 Original code.      #
##################################

import argparse
import gc
import os
import pysam
import sys
import numpy as np
import pandas as pd

__version__ = "version 1.0"

#####################
# Get Error Message #
#####################
def get_error(field, syntax):
    error_msg = "[Error]\nUnknown {} '{}', please check the document for more details.".format(field,syntax)
    return error_msg


#####################################
# Detect Description Field in FASTA #
#####################################
def extract_field(df, name):
    tmp_df = pd.DataFrame(df['read_id'].str.split().to_list())
    target_index = -1
    for i in tmp_df.columns:
        if tmp_df.iloc[0,i][:len(name)] == name:
            target_index = i
            break
    if target_index == -1:
        return tmp_df[0], pd.DataFrame()
    else:
        return tmp_df[0], tmp_df[target_index].str[len(name)+1:]
    

#############################
# Auto-detect and Read Data #
#############################
def read_data(data):
    
    # initialize
    cmt = ""
    error_msg = get_error('format of input', data)

    if data.name=='<stdin>':
        # CSV stdin
        try:
            for line in data:
                if line[0]=='#': cmt += line   
                else: break
            header = line[:-1].split(',')
            df = pd.read_csv(data, comment='#', names=header)
        except:
            print(error_msg)
            sys.exit(1)

    else:
        data = data.name
        filename, file_extension = os.path.splitext(data)
        file_extension = file_extension[1:]

        # CSV
        if file_extension=='csv':
            with open(data, 'r') as f:
                for line in f:
                    if line[0]=='#': cmt += line   
                    else: break
                header = line[:-1].split(',')
                df = pd.read_csv(f, comment='#', names=header)
        
        # TSV
        elif file_extension=='tsv' or file_extension=='bed':
            with open(data, 'r') as f:
                for line in f:
                    if line[0]=='#': cmt += line   
                    else: break
                header = line[:-1].split(',')
                df = pd.read_csv(f, comment='#', sep='\t', names=header)
        
        # FASTA
        elif file_extension=='fasta' or file_extension=='fa':

            i = -1
            id = []
            seq = []
            with open(data, 'r') as f:
                for line in f:
                    if line[0]=='>':
                        id.append(line[1:-1])
                        seq.append('')
                        i += 1
                    else:
                        seq[i] += line[:-1]       
            df = pd.DataFrame(list(zip(id,seq)), columns=['read_id','read_seq'])

            # check if read_count is in read_id
            # must in the format "read_id|read_count"
            df_tmp = df['read_id'].str.split('|', expand=True)
            if 1 in df_tmp.columns:
                df['read_id'] = df_tmp[0]
                df['read_count'] = df_tmp[1].astype(int)
            del df_tmp
            
        # FASTQ
        elif file_extension=='fastq' or file_extension=='fq':

            i = 0
            id = []
            seq = []
            with open(data, 'r') as f:
                for line in f:
                    if i%4==0 and line!='\n':
                        id.append(line[1:-1].split()[0])
                    elif i%4==1 and line!='\n':
                        seq.append(line[:-1])
                    i += 1
            df = pd.DataFrame(zip(id,seq), columns=['read_id','read_seq'])
            df['read_count'] = 1

        # BAM
        elif file_extension=='bam':

            bamfile = pysam.AlignmentFile(data, 'rb')
            fields = {  
                'read_seq': [],    # segment SEQuence
                'read_id': [],     # Query template NAME
                'ref_id': [],      # References sequence NAME
                'target_pos': [],  # 0- based leftmost mapping POSition
                'score': [],       # MAPping Quality
                'CIGAR': [],       # CIGAR string
                'MD': [],          # "MD:Z" in <Bowtie2>
            }

            for read in bamfile.fetch(until_eof=True):
                fields['read_seq'].append(read.query_sequence)
                fields['read_id'].append(read.query_name)
                fields['ref_id'].append(read.reference_name)
                fields['target_pos'].append(read.reference_start)
                fields['score'].append(read.mapping_quality)
                fields['CIGAR'].append(read.cigarstring)
                fields['MD'].append(read.get_tag("MD"))
            bamfile.close()

            df = pd.DataFrame.from_dict(fields)
            
            # get position (0-based)
            df['init_pos'] = df['target_pos'].astype(int) + 1
            df['end_pos'] = df['target_pos'] + df['read_seq'].str.len()
            del df['target_pos']

            # check if read_count is in read_id
            # must in the format "read_id|read_count"
            df_tmp = df['read_id'].str.split('|', expand=True)
            if 1 in df_tmp.columns:
                df['read_id'] = df_tmp[0]
                df['read_count'] = df_tmp[1].astype(float)
            del df_tmp
         
        # SAM
        elif file_extension=='sam':

            sam_dict = {
                # Mandatory Fields
                0: 'read_id',    # Query template NAME
                1: 'FLAG',       # bitwise FLAG
                2: 'ref_id',     # References sequence NAME
                3: 'target_pos', # 1- based leftmost mapping POSition
                4: 'score',      # MAPping Quality
                5: 'CIGAR',      # CIGAR string
                6: 'RNEXT',      # Ref. name of the mate/next read
                7: 'PNEXT',      # Position of the mate/next read
                8: 'TLEN',       # observed Template LENgth
                9: 'read_seq',   # segment SEQuence
                10: 'QUAL',      # ASCII of Phred-scaled base QUALity+33

                # Optional Fiels
                # Below is for <Bowtie2>
                17: 'MD1',
                18: 'MD2',
            }
            main_index = [0,2,3,4,5,9]
            optional_index = [17,18]

            try:
                df = pd.read_csv(data, comment='@', sep='\t', header=None, usecols=main_index+optional_index)
            except:
                df = pd.read_csv(data, comment='@', sep='\t', header=None, usecols=main_index)
            df = df.rename(columns=sam_dict)
            
            # get position (1-based)
            df['init_pos'] = df['target_pos'].astype(int)
            df['end_pos'] = df['target_pos'] + df['read_seq'].str.len() - 1
            del df['target_pos']

            # get MD field
            if 'MD1' in df.columns and 'MD2' in df.columns:
                df.loc[ df['MD1'].str[:2]=='MD', 'MD'] = df['MD1'].str[5:]
                df.loc[ df['MD2'].str[:2]=='MD', 'MD'] = df['MD2'].str[5:]
                del df['MD1']
                del df['MD2']
            
            # check if read_count is in read_id
            # must in the format "read_id|read_count"
            df_tmp = df['read_id'].str.split('|', expand=True)
            if 1 in df_tmp.columns:
                df['read_id'] = df_tmp[0]
                df['read_count'] = df_tmp[1].astype(float)
            del df_tmp

        else:
            print("[Error]")
            print("Unknown file extension {}, which can only be 'fasta/fa', 'fastq/fq', 'tsv/bed', 'csv' or 'sam/bam'.".format(file_extension))
            sys.exit(1)
    
    return df, cmt


##############################################
# Filter Sequence with Length and Nucleotide #
##############################################
def seq_filter(df, syntax):

    # initialize
    cmt = '# filter={}\n'.format(syntax)
    error_msg = get_error('filter', syntax)
    syntax = syntax.split(':')
    seq_len = df['read_seq'].str.len()

    # length
    if syntax[1]=='':
        pass
    elif '-' in syntax[1]:
        lst = syntax[1].split('-')
        num1,num2 = int(lst[0]),int(lst[1])
        df = df[ (seq_len>=num1) & (seq_len<=num2) ].reset_index(drop=True)
    elif '>=' in syntax[1]:
        num = int(syntax[1].split('>=')[1]) 
        df = df[ seq_len >= num ].reset_index(drop=True)
    elif '>' in syntax[1]:
        num = int(syntax[1].split('>')[1])
        df = df[ seq_len > num ].reset_index(drop=True)
    elif '<=' in syntax[1]:
        num = int(syntax[1].split('<=')[1]) 
        df = df[ seq_len <= num ].reset_index(drop=True)
    elif '<' in syntax[1]:
        num = int(syntax[1].split('<')[1])
        df = df[ seq_len < num ].reset_index(drop=True)
    else:
        print(error_msg)
        sys.exit(1)
    
    # nucleotide
    # position = (syntax_pos),(nucleotide_pos) in head and tail
    for position in [[0,0],[2,-1]]:
        if syntax[position[0]]!='':
            TF_table = [False]*len(df)
            for nt in syntax[position[0]]:
                if nt in ['A','T','C','G','U']:
                    TF_table = ((df['read_seq'].str[position[1]]==nt) | TF_table)
                else:
                    print(error_msg)
                    sys.exit(1)
            df = df[ TF_table ].reset_index(drop=True)

    return df, cmt


######################################
# Merge Columns to Read or Reference #
######################################
def merge(df, df_reference, syntax):
    
    # initialize
    cmt = "# merge={}\n".format(syntax)
    error_msg = get_error('syntax', syntax)
    syntax = syntax.split(':')
    syntax[1] = syntax[1].split(',')

    # check which columns will be merged, others will be appended
    merge_col = []
    for column in syntax[1]:
        if column in df.columns:
            merge_col.append(column)
    
    # check if merge to read or reference data, default in 'read'
    if syntax[0]=='read' or syntax[0]=='':
        syntax[1].append('read_id')
        merge_col.append('read_id')
    elif syntax[0]=='ref':
        syntax[1].append('ref_id')
        merge_col.append('ref_id')
        df_reference = df_reference.rename(columns={'read_id':'ref_id'})
    else:
        print(error_msg)
        sys.exit(1)
    
    # reassemble and merge two dataframe
    df_reference = df_reference[syntax[1]]
    df = pd.merge(df, df_reference, on=merge_col)
    
    return df, cmt


########################
# Collapse Read Counts #
########################
def collapse(df):
    
    cmt = "# collapse\n"
    
    if 'read_id' in df.columns:
        del df['read_id']
    df = df.groupby('read_seq').count().reset_index()
    df.index = df.index + 1
    df['read_id'] = 'R' + df.index.astype(str)
    df = df[['read_id','read_seq','read_count']].reset_index(drop=True)
    
    return df, cmt


#############################################
# Distribute Read-Count by Duplicated Reads #
#############################################
def distribute(df):

    cmt = "# distribute\n"

    df_tmp = pd.DataFrame()
    df_tmp['dup'] = df.pivot_table(columns=['read_id'], aggfunc='size')
    df_tmp = df_tmp.reset_index()
    df = pd.merge(df, df_tmp, how='left', on='read_id')
    df['read_count'] = df['read_count']/df['dup']
    del df['dup']

    return df, cmt


######################
# Reverse Complement #
######################
def reverse_complement(df, A_map='T'):
    
    def rc(seq):
        rev_dict = {'A':A_map, A_map:'A','C':'G','G':'C'}
        seq = ''.join(rev_dict[nt] for nt in seq[::-1])
        return seq
    
    cmt = "# reverse_complement\n"
    df['read_seq'] = df['read_seq'].apply(lambda x: rc(x))
    
    return df, cmt


##################################
# Calculate Normalization Factor #
##################################
def norm_factor(df, M=1000000, digit=5):
    factor = round(M/df['read_count'].sum(), digit)
    cmt = "# norm_factor={}\n".format(factor)
    return df, cmt


############################
# Read Count Normalization #
############################
def normalize(df, factor):
    df['read_count'] = df['read_count']*float(factor)
    cmt = "# normalize={}\n".format(factor)
    return df, cmt


##################################
# Output File in Specific Format #
##################################
def output_file(df, output, cmt='', format='csv'):
    
    error_msg = get_error('format', format)
    
    if format=='csv':
        output.write(cmt)
        df.to_csv(output, index=False)

    elif format=='fasta' or format=='fa':
        df['read_id'] = '>' + df['read_id'] + '|' + df['read_count'].astype(str)
        df = df[['read_id','read_seq']]
        df.to_csv(output, index=False, header=None, sep='\n')
    else:
        print(error_msg)
        sys.exit(1)


################
# Main Program #
################
if __name__ == '__main__':

    # arguments
    parser = argparse.ArgumentParser(
        description="This program is a universal tool to preprocess NGS-Seq.",
    )
    parser.add_argument("-v", "--version",
                        action="version",
                        version="%(prog)s " + __version__)
    parser.add_argument("-i", "--input",
                        nargs='?',
                        type=argparse.FileType('r'),
                        default=sys.stdin,
                        help="input file in 'csv', 'tsv', 'fasta', 'fastq' or 'sam' format")
    parser.add_argument("-o", "--output",
                        nargs='?',
                        type=argparse.FileType('w'),
                        default=sys.stdout,
                        help="output file in 'csv' or 'fasta' format")
    parser.add_argument("-r", "--reference",
                        nargs='?',
                        type=argparse.FileType('r'),
                        help="reference data for some operations")
    parser.add_argument("--collapse",
                        action="store_true",
                        help="collapse duplicated reads into one read with accumulated read-count")
    parser.add_argument("--filter",
                        help="filter sequences with length or head's and tail's nucleotide")
    parser.add_argument("--merge",
                        help="merge specific column of reference data to input data")
    parser.add_argument("--distribute",
                        action="store_true",
                        help="distribute read-count by read ID")
    parser.add_argument("--rc",
                        action="store_true",
                        help="reverse complement, A to T in default")
    parser.add_argument("--norm_factor",
                        action="store_true",
                        help="calculate normalization factor in RPM")
    parser.add_argument("--normalize",
                        type=float,
                        help="normalize read-count by a factor")
    parser.add_argument("--format",
                        choices=["fasta", "fa", "csv"],
                        help='transfer file format')
    args = parser.parse_args()

    # main program
    df = pd.DataFrame()
    comment = ""
    
    if args.input!=None:
        df,cmt = read_data(args.input)
        comment += cmt
    else:
        print("[Error]")
        print("Input data is empty, please setting one with stdin or argument '-i'.")
        sys.exit(1) 

    if args.reference!=None:
        df_refernce,cmt = read_data(args.reference)

    if args.collapse:
        df,cmt = collapse(df)
        comment += cmt
    
    if args.filter!=None:
        df,cmt = seq_filter(df, args.filter)
        comment += cmt
    
    if args.merge!=None:
        df,cmt = merge(df, df_refernce, args.merge)
        comment += cmt

    if args.distribute:
        df,cmt = distribute(df)
        comment += cmt

    if args.rc:
        df,cmt = reverse_complement(df, 'T')
        comment += cmt
    
    if args.norm_factor:
        df,cmt = norm_factor(df)
        comment += cmt
    
    if args.normalize!=None:
        df,cmt = normalize(df, args.normalize)
        comment += cmt
    
    if args.format!=None:
        file_extension = args.format
    else:
        file_extension = 'csv'

    output_file(df, args.output, cmt=comment, format=file_extension)
    