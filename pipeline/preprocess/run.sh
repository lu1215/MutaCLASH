#!/bin/bash

## usage example (only using chira algorithm, so regulator and transcript files are not needed)##
# sh run.sh <input file> <data name>
# other trim_galore parameter need to be written in ../../preprocess.conf
## ---------------------------------- ##

## get parameter from command ##
temp_path=$(dirname "$0")/metadata/
input=$1
data_name=$2
echo input = $input
echo metadata_path = $temp_path
## ---------------------------------- ##

## get parameter from preprocess.conf ##
echo "Getting parameter from preprocess.conf"
. ../../preprocess.conf
## ---------------------------------- ##

## trim_galore part ##
echo "Executing trim_galore"
if test "$link" = "" -o "$link" = "None";then
    trim_galore --length ${len} --dont_gzip -o ${temp_path} -q ${trim} --max_length ${slen} ${input} 
else
    trim_galore --length ${len} --dont_gzip -a ${link} -o ${temp_path} -q ${trim} --max_length ${slen} ${input}
fi
## ---------------------------------- ##

## Deduplication ##
echo "Deduplication part"
python Deduplication.py --data_path ${temp_path}/${data_name}_trimmed.fq --data_name $data_name
## ---------------------------------- ##

## output data ##
# if input data format is .fq, no error will dispaly in terminal
mv ${temp_path}"$data_name".fastq_trimming_report.txt output/${data_name}_trimming.log 2>/dev/null || mv ${temp_path}"$data_name".fq_trimming_report.txt output/${data_name}_trimming.log
# delete metadata (.fq file)
rm ${temp_path}"$data_name"_trimmed.fq
## ---------------------------------- ##