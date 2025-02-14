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
# PRG-1 or ALG-1 NGS data in csv format( after processed by MutaCLASH.sh ) in generate_paper_figures/input/
if [ "$input_file" = "PRG-1" ]; then
    input_path="generate_paper_figures/input/PRG-1.csv"
    regulator_file="data/reference/piRNA_WS275.fa"
    transcript_file="data/reference/mRNA_WS275.fa"
    REGION=10/0/-15/-30
    algorithm="pirScan"
    abundance_type="site"
elif [ "$input_file" = "ALG-1" ]; then
    input_path="generate_paper_figures/input/ALG-1.csv"
    regulator_file="data/reference/miRNA_WS275.fa"
    transcript_file="data/reference/mRNA_WS275.fa"
    REGION=200/140/100/60
    algorithm="miRanda"
    abundance_type="abu"
fi

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


# --------------------------

# set environment

cd ..
. ./environment.sh
cd generate_paper_figures

echo "Step1. add abundance"
cd ../pipeline/add_abundance
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
cd ../../generate_paper_figures/

# --------------------------

echo "Step2. collect files"
DIR=${DATA}_$(date +%Y-%m-%d_%H-%M-%S)
mkdir output/${DIR}
cp ../pipeline/${OUTPUT} output/${DIR}/${DATA}.csv
cp ../pipeline/${OUTPUT} output/${DIR}/${DATA}.csv

cp -r ../pipeline/generate_figure/figure output/${DIR}/
cp -r ../pipeline/generate_figure/log output/${DIR}/
cmd_log=output/${DIR}/log/${DATA}_command.log
touch ${cmd_log}
echo Read File: $input_path >> ${cmd_log}
echo Regulator File: $regulator_file >> ${cmd_log}
echo Transcript File: $transcript_file >> ${cmd_log}
echo Tool: $TOOL >> ${cmd_log}
echo Algorithm: $algorithm >> ${cmd_log}
echo Abundance Analysis Type: $abundance_type >> ${cmd_log}

if [ $DEL_META = true ]
then
    rm ../pipeline/preprocess/output/${DATA}*
    rm ../pipeline/chira/${DATA}*
    rm ../pipeline/find_deletion/ALL_output/${DATA}*
    rm ../pipeline/predict_site/scan_output/${DATA}*
    rm ../pipeline/predict_site/mir_output/${DATA}*
    rm ../pipeline/predict_site/up_output/${DATA}*
    rm ../pipeline/data_processing/after_preprocess/${DATA}*
    rm ../pipeline/${OUTPUT}
fi

echo Output: output/${DIR}
echo "Program complete."
