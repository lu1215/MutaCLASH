#!/bin/bash
# ===========================
# Main Arguments
# <input file>: NGS data in FASTQ format.
# <regulator file>: regulator file in FASTA format.
# <transcript file>: transcript file in FASTA format.
# <tool>: Tool used to detect hybrid reads, which can be "chira", "hyb", "clan".
# <algorithm>: Algorithm used to predict binding sites, which can be "pirScan", "miRanda", "RNAup".
# <abundance analysis type>: Method used to analyze abundance, which can be "abu", "region", "site", "up".
# Full Documentation: https://github.com/RyanCCJ/MutaCLASH
# ===========================

# read path
READ=../../$1
# regulator path
REG=../../$2
# target path
TAR=../../$3
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
TOOL=$4
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
    sh run.sh ${DATA} ../preprocess/output/${DATA}.fa ${REG} ${TAR} ${HYBRID} 4 12 6 4 18
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
sh run.sh ${TOOL} ../${BWA_OUTPUT} ../${OUTPUT} ${REG} ${TAR}
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

echo "Step5. data processing"
cd data_processing
# >>>
sh run.sh ../${OUTPUT} ${TAR}
# >>>
OUTPUT=data_processing/after_preprocess/${DATA}_${TOOL}_step1_scan_mir_RNAup_final.csv
cd ..

# --------------------------

echo "Step6. add abundance"
cd add_abundance
# [n/extend_length]
EXTEND=25
# [region/site/up/abu]
if [ -n "$6" ]
then
    TYPE=$6
else
    TYPE=none
fi
# >>>
sh run.sh ../${OUTPUT} ${REG} ${TAR} ${EXTEND} ${TYPE}
# >>>
if [ $TYPE = "abu" ]
then
    OUTPUT=add_abundance/add_abu_info/abu_${EXTEND}_${DATA}_${TOOL}_step1_scan_mir_RNAup_final.csv
elif [ $TYPE = "region" ] || [ $TYPE = "site" ] || [ $TYPE = "up" ]
then
    OUTPUT=add_abundance/add_22g_info/22g_${TYPE}_${EXTEND}_${DATA}_${TOOL}_step1_scan_mir_RNAup_final.csv
fi
cd ..

# --------------------------

echo "Step7. generate figure"
cd generate_figure
# [pirScan/miRanda/RNAup]
Algorithm=$5
# 22G normalization factor
G22_FACTOR=811.03  # WAGO-1_IP WT
# abundance region, leave blank for 2/3 and 1/3
# miRNA: 200/140/100/60
# piRNA: 10/0/-15/-30
# REGION=10/0/-15/-30
# [png/svg]
FIGURE=png
# >>>
sh run.sh ${DATA} ../${OUTPUT} ${Algorithm} ${TYPE} ${G22_FACTOR} ${TAR} ${FIGURE} ${REGION}
# >>>
cd ../../

# --------------------------

echo "Step8. collect files"
DIR=${DATA}_$(date +%Y-%m-%d_%H-%M-%S)
mkdir data/output/${DIR}
cp pipeline/${OUTPUT} data/output/${DIR}/${DATA}.csv
cp -r pipeline/generate_figure/figure data/output/${DIR}/
cp -r pipeline/generate_figure/log data/output/${DIR}/
cp pipeline/clash_analyst/output/${DATA}_trimming.log data/output/${DIR}/log/
cmd_log=data/output/${DIR}/log/${DATA}_command.log
touch ${cmd_log}
echo Read File: $1 >> ${cmd_log}
echo Regulator File: $2 >> ${cmd_log}
echo Transcript File: $3 >> ${cmd_log}
echo Tool: $4 >> ${cmd_log}
echo Algorithm: $5 >> ${cmd_log}
echo Abundance Analysis Type: $6 >> ${cmd_log}

if [ $DEL_META = true ]
then
    rm pipeline/clash_analyst/output/${DATA}*
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
