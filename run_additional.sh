#!/bin/bash
# ===========================
# Main Arguments
# <input file>: ALG-1 or PRG-1.
# Full Documentation: https://github.com/lu1215/MutaCLASH
# ===========================

# example: sh run_additional.sh --input PRG-1
# example: sh run_additional.sh --input ALG-1
# Function to display usage instructions
usage() {
    echo "Usage: $0 --input <input file>"
    exit 1
}

# Parse command-line arguments
while [ $# -gt 0 ]; do
    case $1 in
        --input)
            input_file="$2"
            shift 2
            ;;
        *)
            echo "Unknown parameter: $1"
            usage
            ;;
    esac
done

# Check if required parameters are provided
if [ -z "$input_file" ]; then
    echo "Error: Missing required arguments."
    usage
fi

# Check if input file is valid
if [ "$input_file" != "PRG-1" ] && [ "$input_file" != "ALG-1" ]; then
    echo "Error: Invalid input file. Please provide a valid input file.(PRG-1 or ALG-1)"
    usage
fi

# Default parameters for Tool
TOOL="chira"

# before start the program, please make sure the following files are in the correct path:
# PRG-1 or ALG-1 NGS data in csv format( after processed by MutaCLASH.sh ) in data/input/
if [ "$input_file" = "PRG-1" ]; then
    input_path="data/input/PRG-1.csv"
    regulator_file="data/reference/piRNA_WS275.fa"
    transcript_file="data/reference/mRNA_WS275.fa"
    REGION=10/0/-15/-30
    algorithm="pirScan"
    abundance_type="site"
elif [ "$input_file" = "ALG-1" ]; then
    input_path="data/input/ALG-1.csv"
    regulator_file="data/reference/miRNA_WS275.fa"
    transcript_file="data/reference/mRNA_WS275.fa"
    REGION=200/140/100/60
    algorithm="miRanda"
    abundance_type="abu"
fi

# change back the column names
# sed -i \
#     -e '1s|CLASH read sequence|hybrid_seq|' \
#     -e '1s|read count|read_count|' \
#     -e '1s|Regulator RNA Name|regulator_name|' \
#     -e '1s|Target RNA Name|transcript_name|' \
#     -e '1s|Target RNA Region Found in CLASH Read|rem_tran_target_pos|' \
#     -e '1s|Regulator RNA Region Found in CLASH Read|reg_hyb_target_pos|' \
#     -e '1s|Region on CLASH Read identified as Regulator RNA|on_reg_pos|' \
#     -e '1s|Region on CLASH Read identified as Target RNA|remain_pos|' \
#     -e '1s|pirScan score|targeting_score|' \
#     -e '1s|miRanda score|mir_score|' \
#     -e '1s|Extended Clash Identified Region Start Position (miRanda)|mir_init_pos|' \
#     -e '1s|Extended Clash Identified Region End Position (miRanda)|mir_end_pos|' \
#     -e '1s|miRanda Defined Binding Region (Relative to Extended Clash Identified Region)|mir_target_pos|' \
#     -e '1s|Transcript Binding Sequence (miRanda)|mir_transcript_seq|' \
#     -e '1s|Regulator Binding Sequence (miRanda)|mir_regulator_seq|' \
#     -e '1s|Extended Clash Identified Region Start Position (RNAup)|up_init_pos|' \
#     -e '1s|Extended Clash Identified Region End Position (RNAup)|up_end_pos|' \
#     -e '1s|Transcript Binding Sequence (RNAup)|RNAup_transcript_seq|' \
#     -e '1s|Regulator Binding Sequence (RNAup)|RNAup_regulator_seq|' \
#     -e '1s|RNAup Defined Binding Region (Relative to Clash Identified Region)|RNAup_target_pos|' \
#     -e '1s|RNAup Binding Energy|RNAup_score|' \
#     -e '1s|Deletion Sites on mRNA (Absolute Positions)|D|' \
#     -e '1s|Mismatch Sites on mRNA (Absolute Positions)|M|' \
#     -e '1s|Site-Level Preprocessing (Read Count = 1)|count|' \
#     -e '1s|Normalized Read Count (After Read Deduplication)|Nor_readcount|' \
#     -e '1s|Normalized Count (After Read Deduplication)|Nor_count|' \
#     -e '1s|Overlapping Region Between Regulator and Transcript (Hybrid Read Coordinates)|Overlap|' \
#     -e '1s|mRNA Length|mRNA_len|' \
#     -e '1s|Transcript-Regulator Pair (For Pair Counting)|Hybrid_read|' \
#     -e '1s|,Mutation Sites on mRNA (Deletion + Mismatch, Absolute Positions)|,A|' \
#     $input_path

# read path
READ=../../$input_path
# regulator path
REG=../../$regulator_file
# target path
TAR=../../$transcript_file
# data base name
DATA=$(basename ${READ})
DATA=${DATA%.*}
REG=${REG%.*}.csv
TAR=${TAR%.*}.csv

# remove metadatas
DEL_META=false

# set environment
. ./environment.sh

# --------------------------

echo "Step1. add abundance"
cd pipeline/add_abundance
# [n/extend_length]
EXTEND=25
# [region/site/up/abu]
if [ -n "$abundance_type" ]
then
    TYPE=$abundance_type
else
    TYPE=none
fi
# >>>
sh run.sh ../../${input_path} ${REG} ${TAR} ${EXTEND} ${TYPE}
# >>>
if [ $TYPE = "abu" ]
then
    OUTPUT=add_abundance/add_abu_info/abu_${EXTEND}_${DATA}.csv
elif [ $TYPE = "region" ] || [ $TYPE = "site" ] || [ $TYPE = "up" ]
then
    OUTPUT=add_abundance/add_22g_info/22g_${TYPE}_${EXTEND}_${DATA}.csv
fi
cd ..

# --------------------------

echo "Step2. generate figure"
cd generate_figure
# [pirScan/miRanda/RNAup]
Algorithm=$algorithm
# 22G normalization factor
G22_FACTOR=811.03  # WAGO-1_IP WT
# abundance region, leave blank for 2/3 and 1/3
# miRNA: 200/140/100/60
# piRNA: 10/0/-15/-30
# REGION=10/0/-15/-30
# [png/svg]
FIGURE=svg
# >>>
sh run.sh ${DATA} ../${OUTPUT} ${Algorithm} ${TYPE} ${G22_FACTOR} ${TAR} ${FIGURE} ${REGION}
# >>>
cd ../../

# --------------------------

echo "Step2. collect files"
DIR=${DATA}_$(date +%Y-%m-%d_%H-%M-%S)
mkdir data/output/${DIR}
cp pipeline/${OUTPUT} data/output/${DIR}/${DATA}.csv
cp pipeline/${OUTPUT} data/output/${DIR}/${DATA}.csv
# rename name of columns in data/output/${DIR}/${DATA}.csv
cp data/output/${DIR}/${DATA}.csv data/output/${DIR}/${DATA}_short.csv

python FilterReorderCsv.py data/output/${DIR}/${DATA}_short.csv data/output/${DIR}/${DATA}_short.csv

# rename name of columns in data/output/${DIR}/${DATA}.csv
sed -i \
    -e '1s|hybrid_seq|CLASH read sequence|' \
    -e '1s|read_count|read count|' \
    -e '1s|regulator_name|Regulator RNA Name|' \
    -e '1s|transcript_name|Target RNA Name|' \
    -e '1s|rem_tran_target_pos|Target RNA Region Found in CLASH Read|' \
    -e '1s|reg_hyb_target_pos|Regulator RNA Region Found in CLASH Read|' \
    -e '1s|on_reg_pos|Region on CLASH Read identified as Regulator RNA|' \
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
    data/output/${DIR}/${DATA}_short.csv

cp -r pipeline/generate_figure/figure data/output/${DIR}/
cp -r pipeline/generate_figure/log data/output/${DIR}/
# cp pipeline/preprocess/output/${DATA}_trimming.log data/output/${DIR}/log/
cmd_log=data/output/${DIR}/log/${DATA}_command.log
touch ${cmd_log}
echo Read File: $input_path >> ${cmd_log}
echo Regulator File: $regulator_file >> ${cmd_log}
echo Transcript File: $transcript_file >> ${cmd_log}
echo Tool: $TOOL >> ${cmd_log}
echo Algorithm: $algorithm >> ${cmd_log}
echo Abundance Analysis Type: $abundance_type >> ${cmd_log}

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
