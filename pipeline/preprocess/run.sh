#!/bin/bash

## usage example (only using chira algorithm, so regulator and transcript files are not needed)##
## sh run.sh <input file> <data name> <prep_tool> <dedup>
## <prep_tool>: "cutadapt" or "trim_galore" (cutadapt or others)
## <dedup>: "hyb" or "python" (hyb or others)
## other trim_galore parameter need to be written in ../../preprocess.conf
## ---------------------------------- ##

## get parameter from command ##
temp_path=$(dirname "$0")/metadata/
input=$1
data_name=$2
prep_tool=$3
dedup=$4
echo input = $input
echo metadata_path = $temp_path

## ---------------------------------- ##

## get parameter from preprocess.conf ##
echo "Getting parameter from preprocess.conf"
. ../../preprocess.conf
## ---------------------------------- ##

## preprocessing ##
if [ "$prep_tool" = "cutadapt" ]; then
    ## -------------------- cutadapt -------------------- ##
    echo "Executing cutadapt"
    if test "$link" = "" -o "$link" = "None";then
        echo "without adapter information, stop executing" >&2
        exit 1
    else
        report="${temp_path%/}/${data_name}.fastq_trimming_report.txt"
        cutadapt -m ${len} -M ${slen} -q ${trim} -a ${link} -o ${temp_path}/${data_name}_trimmed.fq ${input} >"$report" 2>&1
    fi
    ## ------------------------------------------------- ##
else
    ## ------------------- trim_galore ------------------ ##
    echo "Executing trim_galore"
    if test "$link" = "" -o "$link" = "None";then
       trim_galore --length ${len} --dont_gzip -o ${temp_path} -q ${trim} --max_length ${slen} ${input}
    else
        trim_galore --length ${len} --dont_gzip -a ${link} -o ${temp_path} -q ${trim} --max_length ${slen} ${input}
    fi
    ## ------------------------------------------------- ##
fi


## trim_galore part ##
# echo "Executing trim_galore"
# if test "$link" = "" -o "$link" = "None";then
#     trim_galore --length ${len} --dont_gzip -o ${temp_path} -q ${trim} --max_length ${slen} ${input} 
# else
#     trim_galore --length ${len} --dont_gzip -a ${link} -o ${temp_path} -q ${trim} --max_length ${slen} ${input}
# fi
## ---------------------------------- ##

## cutadapt ##
# echo "Executing cutadapt"
# if test "$link" = "" -o "$link" = "None";then
#     echo "without adapter information, stop executing" >&2
#     exit 1
# else
#     report="${temp_path%/}/${data_name}.fastq_trimming_report.txt"
#     cutadapt -m ${len} -M ${slen} -q ${trim} -a ${link} -o ${temp_path}/${data_name}_trimmed.fq ${input} >"$report" 2>&1
# fi
## ---------------------------------- ##

## deduplication part ##
if [ "$dedup" = "hyb" ]; then
    echo "Deduplication part (hyb)"
    ../../bio_tool/clash_analyst/hyb-master/bin/solexa2fasta.awk ${temp_path}/${data_name}_trimmed.fq \
    | ../../bio_tool/clash_analyst/hyb-master/bin/fasta2tab.awk \
    > "${temp_path%/}/${data_name}_trimmed.tab"

    ../../bio_tool/clash_analyst/hyb-master/bin/make_comp_fasta.pl \
    "${temp_path%/}/${data_name}_trimmed.tab" \
    > "output/${data_name}.fa"

    rm ${temp_path%/}/${data_name}_trimmed.tab
else
    echo "Deduplication part (python)"
    python Deduplication.py --data_path ${temp_path}/${data_name}_trimmed.fq --data_name "$data_name"
fi

# ## Deduplication ##
# echo "Deduplication part (python)"
# python Deduplication.py --data_path ${temp_path}/${data_name}_trimmed.fq --data_name $data_name

# echo "Deduplication part (hyb)"
# ../../bio_tool/clash_analyst/hyb-master/bin/solexa2fasta.awk ${temp_path}/${data_name}_trimmed.fq | ../../bio_tool/clash_analyst/hyb-master/bin/fasta2tab.awk > ${temp_path}/${data_name}_trimmed.tab
# ../../bio_tool/clash_analyst/hyb-master/bin/make_comp_fasta.pl ${temp_path}/${data_name}_trimmed.tab > ${temp_path}/${data_name}.fasta
# ## ---------------------------------- ##

## output data ##
# if input data format is .fq, no error will dispaly in terminal
mv ${temp_path}"$data_name".fastq_trimming_report.txt output/${data_name}_trimming.log 2>/dev/null || mv ${temp_path}"$data_name".fq_trimming_report.txt output/${data_name}_trimming.log
# delete metadata (.fq file)
rm ${temp_path}"$data_name"_trimmed.fq
## ---------------------------------- ##