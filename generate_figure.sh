# ex: sh add_abundance.sh <read name> <input file> <transcript file> <algorithm> <regulator type> <abundance analysis type> <figure type>

# read name
DATA=$1
# input file base name
INDATA=$(basename ${2})
INDATA=${INDATA%.csv}
# target path
TAR=../../$4
TAR=${TAR%.*}.csv

# remove metadatas
DEL_META=false

# set environment
. ./environment.sh

echo "Step1. generate figure"
cd generate_figure
# [pirScan/miRanda/RNAup]
Algorithm=$5
# 22G normalization factor
G22_FACTOR=811.03  # WAGO-1_IP WT
# abundance region, leave blank for 2/3 and 1/3
# miRNA: 200/140/100/60
# piRNA: 10/0/-15/-30
miRNA_region="200/140/100/60"
piRNA_region="10/0/-15/-30"

if [[ "$6" == "miRNA" ]]; then
    REGION="200/140/100/60"
elif [[ "$6" == "piRNA" ]]; then
    REGION="10/0/-15/-30"
else
    REGION=""
fi

TYPE=$6

# REGION=10/0/-15/-30
# [png/svg]
if [[ "$7" == "svg" ]]; then
    FIGURE="svg"
else
    FIGURE="png"
fi

cd 
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
echo Read File: $1 >> ${cmd_log}
echo Regulator File: $2 >> ${cmd_log}
echo Transcript File: $3 >> ${cmd_log}
echo Tool: $4 >> ${cmd_log}
echo Algorithm: $5 >> ${cmd_log}
echo Abundance Analysis Type: $6 >> ${cmd_log}

echo Output: data/output/${DIR}
echo "Program complete."