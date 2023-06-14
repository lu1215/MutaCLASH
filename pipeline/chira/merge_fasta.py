import os
import pandas as pd 

# parameters
base = os.getcwd() + '/input_data/'

INPUT_DATA = [base + 'HOMO_CLIP1.fa',
	      base + 'HOMO_CLIP2.fa',
              base + 'HOMO_CLIP3.fa',
              base + 'HOMO_CLIP4.fa',
            ]

OUTPUT_DATA = base + 'HOMO_CLIP1-4.fa'


# functions
def fasta_to_csv(data, sep=' ', columns=None):
    name,seq = [],[] 
    with open(data, 'r') as f:
        lines = f.read().splitlines()
        for line in lines:
            if line[0]=='>':
                name.append(line[1:])
            else:
                seq.append(line)
    df1 = pd.DataFrame(zip(name,seq), columns =['name','sequence'])
    df2 = df1['name'].str.split(sep, expand=True)
    df = pd.concat([df1, df2], axis=1)
    
    del df['name']
    del df1
    del df2
    
    if columns!=None:
        for i in range(len(columns)):
            df.rename(columns={i: columns[i]}, inplace=True)
  
    return df


def csv_to_fasta(data, output):
    
    if type(data)==str:
        df = pd.read_csv(data) 
    else:
        df = data.copy()
    
    df = df[['sequence', 'read_count']] 
    df['index'] = df.index.astype(str)
    df['name'] = '>' + df['index'] + '_' + df['read_count'].astype(str)
    df = df[['name','sequence']]
    
    df.to_csv(output, index=False, header=False, sep='\n')


# main
if __name__ == '__main__':

    print("Start reading datas...")
    df_lst = []
    for data in INPUT_DATA:
        print(data)
        df = fasta_to_csv(data, sep='_', columns=['name','read_count'])
        df_lst.append(df)

    print("Start merging datas...")
    df = pd.concat(df_lst, ignore_index=True)
    del df['name']
    df['read_count'] = df['read_count'].astype(int)
    df = df.groupby('sequence')[['read_count']].sum().reset_index()

    print("Start outputting datas...")
    print(OUTPUT_DATA)
    csv_to_fasta(df, OUTPUT_DATA)

    print("Program end with success.")

