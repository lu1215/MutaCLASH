#!/bin/bash

shell_folder=$(cd "$(dirname "$0")";pwd)
. ${shell_folder}/framefunction.sh

# transfer to .fasta
# ---------------------------------------
function csvTofasta(){
    echo "transfer $1 to fasta"
    local filename=$1
    local id_name=$2
    local outputfile=`echo $filename | sed "s/\.csv$/temp\.fasta/"`

    # find sequence column(from framefunction.sh)
    find_col $3 $filename
    local sequence_col=$?

    # beacause first row is column name,we need to ignore it
    # id_name need to use NR-1
    awk -F, -v col=$sequence_col -v id=$id_name 'NR==1 {next}{printf ">%s%s\n%s\n",id,NR-1,$col}' $filename | sed "s/U/T/g" > $outputfile

}
#csvTofasta ../RNA/biotype=mRNA.csv RNA sequence 


# create a id,seq,seq_len file 
# replace sequence U to T
# ---------------------------------------
function addID(){
    echo "add id and length to $1" 
    local filename=$1
    local id_name=$2
    local outputfile=$3

    # find sequence column(from framefunction.sh)
    find_col sequence $filename
    local sequence_col=$?


    head -n1 $filename|grep -q "${id_name}0"
    fd=$?
    if [ "$fd" = "1" ];then
        # column : id,col_seq,length
        awk -F, -v id=$id_name -v col=$sequence_col '{if (NR==1){sub($col,id"_seq",$col);print id""(NR-1),$0,id"_len";} else {print id""(NR-1),$0,length($col);}}' OFS="," $filename | sed "s/U/T/g" > $outputfile
    else
        cp $filename $outputfile
        sed -i "1s/sequence/${id_name}_seq/1" $outputfile
    fi
}

