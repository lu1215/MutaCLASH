#!/bin/bash
# ===========================
# Main Arguments
# <read name>: the name of NGS data
# <input file>: NGS data after processed by MutaCLASH.sh and add_abundance.sh in CSV format.
# <transcript file>: transcript file in CSV format.
# <algorithm>: Algorithm used to predict binding sites, which can be "pirScan", "miRanda", "RNAup".
# <abundance analysis type>: Method used to analyze abundance, which can be "abu", "region", "site", "up".
# Full Documentation: https://github.com/RyanCCJ/MutaCLASH
# ===========================
# ex: sh generate_figure.sh <read name> <input file> <transcript file> <algorithm> <abundance analysis type>
# sh generate_figure.sh SRR6512653.1 data/output/SRR6512653.1_2025-01-26_16-47-30/abu_25_SRR6512653.1.csv data/reference/mRNA_WS275.fa miRanda abu

# read name
DATA=$1
# input file base name
INDATA=$(basename ${2})
INDATA=${INDATA%.csv}
# target path
TAR=../../$3
TAR=${TAR%.*}.csv

# remove metadatas
DEL_META=false

# set environment
. ./environment.sh

echo "Step1. generate figure"
cd pipeline/generate_figure
# [pirScan/miRanda/RNAup]
Algorithm=$4
# 22G normalization factor
G22_FACTOR=811.03  # WAGO-1_IP WT
# abundance region, leave blank for 2/3 and 1/3
# abundance region, leave blank for 2/3 and 1/3
# miRNA: 200/140/100/60
# piRNA: 10/0/-15/-30
REGION=200/140/100/60

TYPE=$5

# REGION=10/0/-15/-30
# [png/svg]
FIGURE=png

# >>>
sh run.sh ${DATA} ../../$2 ${Algorithm} ${TYPE} ${G22_FACTOR} ${TAR} ${FIGURE} ${REGION}
# >>>
cd ../../

echo "Step2. collect files"
DIR=${DATA}_$(date +%Y-%m-%d_%H-%M-%S)_figure
mkdir data/output/${DIR}
cp -r pipeline/generate_figure/figure data/output/${DIR}/
cp -r pipeline/generate_figure/log data/output/${DIR}/
cmd_log=data/output/${DIR}/log/${DATA}_command.log
touch ${cmd_log}
echo Read File: $2 >> ${cmd_log}
echo Transcript File: $3 >> ${cmd_log}
echo Algorithm: $4 >> ${cmd_log}
echo Abundance Analysis Type: $5 >> ${cmd_log}

echo Output: data/output/${DIR}
echo "Program complete."