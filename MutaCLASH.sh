#!/bin/bash
# ===========================
# Main Arguments
# <input file>: NGS data in FASTQ format.
# <regulator file>: regulator file in FASTA format.
# <transcript file>: transcript file in FASTA format.
# optional arguments:
# <len>: Minimum hybrid length (default: 17).
# <slen>: Maximum hybrid length (default: 70).
# <link>: Adapter sequence (default: "None").
# <trim>: Phred score (default: 30).
# Full Documentation: https://github.com/RyanCCJ/MutaCLASH
# ===========================

# Default preprocessing parameter values
len=17
slen=70
link="None"
trim=30

# default non-transposon mode
transposon=false

# Function to display usage instructions
usage() {
    echo "Usage: $0 --input <input file> --regulator <regulator file> --transcript <transcript file> [--len <min hybrid length>] [--slen <max hybrid length>] [--link <adapter sequence>] [--trim <phred score>] [--transposon]"
    exit 1
}

# Parse command-line arguments
while [ $# -gt 0 ]; do
    case $1 in
        --input)
            input_file="$2"
            shift 2
            ;;
        --regulator)
            regulator_file="$2"
            shift 2
            ;;
        --transcript)
            transcript_file="$2"
            shift 2
            ;;
        --len)
            len="$2"
            shift 2
            ;;
        --slen)
            slen="$2"
            shift 2
            ;;
        --link)
            link="$2"
            shift 2
            ;;
        --trim)
            trim="$2"
            shift 2
            ;;
        --transposon)  
            transposon=true;    
            shift 1 
            ;;
        *)
            echo "Unknown parameter: $1"
            usage
            ;;
    esac
done

# write preprocessing parameter to preprocess.conf
# Use sed to update variables in-place, or append if not found
sed -i "/^len=/c\len=$len" preprocess.conf || echo "len=$len" >> preprocess.conf
sed -i "/^slen=/c\slen=$slen" preprocess.conf || echo "slen=$slen" >> preprocess.conf
sed -i "/^link=/c\link=$link" preprocess.conf || echo "link=$link" >> preprocess.conf
sed -i "/^trim=/c\trim=$trim" preprocess.conf || echo "trim=$trim" >> preprocess.conf

# Check if required parameters are provided
if [ -z "$input_file" ] || [ -z "$regulator_file" ] || [ -z "$transcript_file" ]; then
    echo "Error: Missing required arguments."
    usage
fi

# read path
READ=../../$input_file
# regulator path
REG=../../$regulator_file
# target path
TAR=../../$transcript_file
# data base name
DATA=$(basename ${READ})
DATA=${DATA%.*}

# remove metadatas
DEL_META=false

# set environment
. ./environment.sh

# ===========================

# echo "Step1. clash analyst"
# cd pipeline/clash_analyst
# # [hyb/clan/chira]
# TOOL=$4
# # >>>
# sh run.sh ${READ} ${REG} ${TAR} ${TOOL} ${DATA}
# # >>>
# cd ..
# OUTPUT=clash_analyst/output/${DATA}_${TOOL}.csv

# # --------------------------

echo "Step1. Preprocess(Trim_galore and De-duplication)"
cd pipeline/preprocess
# [hyb/clan/chira]
TOOL="chira"
# >>>
sh run.sh ${READ} ${DATA}
# >>>
cd ..

# --------------------------

echo "Step2. chira"

# [single/chimeras]
HYBRID=chimeras

if [ $TOOL = "chira" ]
then
    cd chira
    # run.sh [data_name] [read] [regulator] [target] [hybrid(chimeras)] [thread(4)] [seed_length(12)] [gap_penalty(6)] [mismatch_penalty(4)] [score_cutoff(18)]
    # >>>
    sh run.sh ${DATA} ../preprocess/output/${DATA}.fa ${REG} ${TAR} ${HYBRID} 4 12 6 4 18 ${transposon}
    # >>>
    cd ..
    TOOL=chira_${HYBRID}
    BWA_OUTPUT=chira/${DATA}_map_dir/sorted.bam
    OUTPUT=chira/${DATA}_extract_dir/${DATA}_${TOOL}.csv
fi

# --------------------------

echo "Step3. find deletion"
cd find_deletion
# >>>
sh run.sh ${TOOL} ../${BWA_OUTPUT} ../${OUTPUT} ${REG} ${TAR} ${transposon}
# >>>
cd ..
OUTPUT=find_deletion/ALL_output/${DATA}_${TOOL}_step1.csv

# --------------------------

echo "Step4. predict site"
cd predict_site
REG=${REG%.*}.csv
TAR=${TAR%.*}.csv
# [n/extend_length]
EXTEND=n

# pirScan
# >>>
sh run_pirScan.sh ../${OUTPUT} ${REG} ${TAR} ${EXTEND}
# >>>
OUTPUT=predict_site/scan_output/${DATA}_${TOOL}_step1_scan.csv

# miRanda
# >>>
sh run_miRanda.sh ../${OUTPUT} ${REG} ${TAR} ${EXTEND}
# >>>
OUTPUT=predict_site/mir_output/${DATA}_${TOOL}_step1_scan_mir.csv

# RNAup
# >>>
sh run_RNAup.sh ../${OUTPUT} ${REG} ${TAR} ${EXTEND}
# >>>
OUTPUT=predict_site/up_output/${DATA}_${TOOL}_step1_scan_mir_RNAup.csv
cd ..

# --------------------------

echo "Step5. single read processing & analysing"
cd single_read_stat

# >>>
sh run_single_read_stat.sh ${DATA} ../../$regulator_file ../../$transcript_file ${OUTPUT}
# >>>
OUTPUT=single_read_stat/output/${DATA}_${TOOL}_step1_scan_mir_RNAup_with_sgl_mut_stat.csv
cd ..
# --------------------------

echo "Step6. data processing"
cd data_processing
# >>>
sh run.sh ../${OUTPUT} ${TAR}
# >>>
OUTPUT=data_processing/after_preprocess/${DATA}_${TOOL}_step1_scan_mir_RNAup_with_sgl_mut_stat_final.csv
cd ..

# --------------------------



echo "Step7. collect files"
cd ..
DIR=${DATA}_$(date +%Y-%m-%d_%H-%M-%S)_MutaCLASH
mkdir data/output/${DIR}
mkdir data/output/${DIR}/log
cp pipeline/${OUTPUT} data/output/${DIR}/${DATA}.csv
cp data/output/${DIR}/${DATA}.csv data/output/${DIR}/${DATA}_short.csv

python FilterReorderCsv.py data/output/${DIR}/${DATA}_short.csv data/output/${DIR}/${DATA}_short.csv

# rename name of columns in data/output/${DIR}/${DATA}.csv
sed -i \
    -e '1s|hybrid_seq|CLASH read sequence|' \
    -e '1s|read_count|read count|' \
    -e '1s|nor_readcount|normalized read count (evenly distributed)|' \
    -e '1s|regulator_name|Regulator RNA Name|' \
    -e '1s|transcript_name|Target RNA Name|' \
    -e '1s|rem_tran_target_pos|Target RNA Region Found in CLASH Read|' \
    -e '1s|reg_hyb_target_pos|Region on CLASH Read identified as Regulator RNAd|' \
    -e '1s|on_reg_pos|Regulator RNA Region Found in CLASH Read|' \
    -e '1s|remain_pos|Region on CLASH Read identified as Target RNA|' \
    -e '1s|targeting_score|pirScan score|' \
    -e '1s|mir_score|miRanda score|' \
    -e '1s|mir_init_pos|Extended Clash Identified Region Start Position (miRanda)|' \
    -e '1s|mir_end_pos|Extended Clash Identified Region End Position (miRanda)|' \
    -e '1s|mir_target_pos|miRanda Defined Binding Region (Relative to Extended Clash Identified Region)|' \
    -e '1s|mir_transcript_seq|Target Binding Sequence (miRanda)|' \
    -e '1s|mir_regulator_seq|Regulator Binding Sequence (miRanda)|' \
    -e '1s|up_init_pos|Extended Clash Identified Region Start Position (RNAup)|' \
    -e '1s|up_end_pos|Extended Clash Identified Region End Position (RNAup)|' \
    -e '1s|RNAup_transcript_seq|Target Binding Sequence (RNAup)|' \
    -e '1s|RNAup_regulator_seq|Regulator Binding Sequence (RNAup)|' \
    -e '1s|RNAup_target_pos|RNAup Defined Binding Region (Relative to Clash Identified Region)|' \
    -e '1s|RNAup_score|RNAup Binding Energy|' \
    -e '1s|D|Deletion Sites on mRNA (Absolute Positions)|' \
    -e '1s|M|Mismatch Sites on mRNA (Absolute Positions)|' \
    -e '1s|Nor_readcount|Normalized Read Count (After Read Deduplication)|' \
    -e '1s|Nor_count|Normalized Count (After Read Deduplication)|' \
    -e '1s|Overlap|Overlapping Region Between Regulator and Transcript (Hybrid Read Coordinates)|' \
    -e '1s|mRNA_len|mRNA Length|' \
    -e '1s|Hybrid_read|Transcript-Regulator Pair (For Pair Counting)|' \
    -e '1s|,A|,Mutation Sites on mRNA (Deletion + Mismatch, Absolute Positions)|' \
    -e '1s|mir_energy|Binding Energy Calculated by miRanda|' \
    -e '1s|pirscan|pirScan|' \
    -e '1s|miranda|miRanda|' \
    -e '1s|rnaup|RNAup|' \
    -e '1s|raw_regulator_seq|pirscan Regulator RNA sequence|' \
    data/output/${DIR}/${DATA}_short.csv

cp pipeline/preprocess/output/${DATA}_trimming.log data/output/${DIR}/log/
cmd_log=data/output/${DIR}/log/${DATA}_command.log
touch ${cmd_log}
echo Read File: $input_file >> ${cmd_log}
echo Regulator File: $regulator_file >> ${cmd_log}
echo Transcript File: $transcript_file >> ${cmd_log}
echo Tool: $TOOL >> ${cmd_log}

if [ $DEL_META = true ]
then
    rm pipeline/preprocess/output/${DATA}*
    rm pipeline/chira/${DATA}*
    rm pipeline/find_deletion/ALL_output/${DATA}*
    rm pipeline/predict_site/scan_output/${DATA}*
    rm pipeline/predict_site/mir_output/${DATA}*
    rm pipeline/predict_site/up_output/${DATA}*
    rm pipeline/data_processing/after_preprocess/${DATA}*
    rm pipeline/${OUTPUT}
fi

echo Output: data/output/${DIR}
echo "Program complete."