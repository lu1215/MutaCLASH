#!/bin/bash
# ===========================
# Main Arguments
# <read name>: the name of NGS data
# <input file>: NGS data after processed by MutaCLASH.sh in csv format.
# <regulator file>: regulator file in CSV format.
# <transcript file>: transcript file in CSV format.
# <abundance analysis type>: Method used to analyze abundance, which can be "abu", "region", "site", "up".
# Full Documentation: https://github.com/RyanCCJ/MutaCLASH
# ===========================
# ex: sh add_abundance.sh <read name> <input file> <regulator file> <transcript file> <abundance analysis type>

# read name
DATA=$1
# regulator path
REG=../../$3
# target path
TAR=../../$4

REG=${REG%.*}.csv
TAR=${TAR%.*}.csv

# input file base name
INDATA=$(basename ${2})
echo $INDATA
INDATA=${INDATA%.csv}
echo $INDATA

# remove metadatas
DEL_META=false

# set environment
. ./environment.sh

echo "Step1. add abundance"
cd pipeline/add_abundance

# [n/extend_length]
EXTEND=25

# [region/site/up/abu]
if [ -n "$5" ]
then
    TYPE=$5
else
    TYPE=none
fi

# >>>
sh run.sh ../../$2 ${REG} ${TAR} ${EXTEND} ${TYPE}
# >>>

if [ $TYPE = "abu" ]
then
    OUTPUT=add_abundance/add_abu_info/abu_${EXTEND}_${INDATA}.csv
elif [ $TYPE = "region" ] || [ $TYPE = "site" ] || [ $TYPE = "up" ]
then
    OUTPUT=add_abundance/add_22g_info/22g_${TYPE}_${EXTEND}_${INDATA}.csv
fi
cd ..

echo "Step2. collect files"
cd ..
echo $(pwd)
DIR=${DATA}_$(date +%Y-%m-%d_%H-%M-%S)_abundance
mkdir data/output/${DIR}
mkdir data/output/${DIR}/log
cp pipeline/${OUTPUT} data/output/${DIR}/
cmd_log=data/output/${DIR}/log/${DATA}_command.log
touch ${cmd_log}
echo Read File: $1 >> ${cmd_log}
echo Regulator File: $2 >> ${cmd_log}
echo Transcript File: $3 >> ${cmd_log}
echo Abundance Analysis Type: $6 >> ${cmd_log}

if [ $DEL_META = true ]
then
    rm pipeline/${OUTPUT}
fi

echo Output: data/output/${DIR}
echo "Program complete."