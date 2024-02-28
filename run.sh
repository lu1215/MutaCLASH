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

# single or chimeras
HYBRID=chimeras

# ===========================

echo "Step1. data checking"
cd pipeline/data_checking
# >>>
sh run.sh ${READ} ${REG} ${TAR}
# >>>
cd ..

# --------------------------

echo "Step2. chira"
cd chira
READ=${READ%.*}.fa
REG=${REG%.*}.fa
TAR=${TAR%.*}.fa
# run.sh [data_name] [read] [regulator] [target] [hybrid(chimeras)] [thread(8)] [seed_length(12)] [gap_penalty(6)] [mismatch_penalty(4)] [score_cutoff(18)]
# >>>
sh run.sh ${DATA} ${READ} ${REG} ${TAR} ${HYBRID} 8 12 6 4 18
# >>>
cd ..
BWA_OUTPUT=chira/${DATA}_map_dir/sorted.bam
OUTPUT=chira/${DATA}_extract_dir/${DATA}_${HYBRID}.csv

# --------------------------

echo "Step3. find deletion"
cd find_deletion

# Bowtie2 (v1)
# TAR=${TAR%.*}.csv
# sh run.sh ../${OUTPUT} ${TAR}

# Samtool (v2) > recommend!
# >>>
sh run_v2.sh ../${BWA_OUTPUT} ../${OUTPUT} ${REG} ${TAR}
# >>>
cd ..
OUTPUT=find_deletion/ALL_output/${DATA}_step1.csv

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
OUTPUT=predict_site/scan_output/${DATA}_step1_scan.csv

# miRanda
# >>>
sh run_miRanda.sh ../${OUTPUT} ${REG} ${TAR} ${EXTEND}
# >>>
OUTPUT=predict_site/mir_output/${DATA}_step1_scan_mir.csv

# RNAup
# >>>
sh run_RNAup.sh ../${OUTPUT} ${REG} ${TAR} ${EXTEND}
# >>>
OUTPUT=predict_site/up_output/${DATA}_step1_scan_mir_RNAup.csv
cd ..

# --------------------------

echo "Step5. data processing"
cd data_processing
# >>>
sh run.sh ../${OUTPUT} ${TAR}
# >>>
OUTPUT=data_processing/after_preprocess/${DATA}_step1_scan_mir_RNAup_final.csv
cd ..

# --------------------------

echo "Step6. abundance"
cd add_abundance
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
sh run.sh ../${OUTPUT} ${REG} ${TAR} ${EXTEND} ${TYPE}
# >>>
if [ $TYPE = "abu" ]
then
    OUTPUT=add_abundance/add_abu_info/abu_${EXTEND}_${DATA}_step1_scan_mir_RNAup_final.csv
elif [ $TYPE = "region" ] || [ $TYPE = "site" ] || [ $TYPE = "up" ]
then
    OUTPUT=add_abundance/add_22g_info/22g_${TYPE}_${EXTEND}_${DATA}_step1_scan_mir_RNAup_final.csv
fi
cd ..

# --------------------------

echo "Step7. generate figure"
cd generate_figure
# [pirScan/miRanda/RNAup]
TOOL=$4
# normalization factor
G22_FACTOR=811.03  # WAGO-1_IP WT
# >>>
sh run.sh ${DATA} ../${OUTPUT} ${TOOL} ${TYPE} ${G22_FACTOR} ${TAR}
# >>>
cd ../../

# --------------------------

echo "Step8. collect files"
DIR=${DATA}_$(date +%Y-%m-%d_%H-%M-%S)
mkdir data/output/${DIR}
cp pipeline/${OUTPUT} data/output/${DIR}/${DATA}.csv
cp -r pipeline/generate_figure/figure data/output/${DIR}/
cp -r pipeline/generate_figure/log data/output/${DIR}/
cmd_log=data/output/${DIR}/log/${DATA}_command.log
touch ${cmd_log}
echo Read File: $1 >> ${cmd_log}
echo Regulator File: $2 >> ${cmd_log}
echo Transcript File: $3 >> ${cmd_log}
echo Algorithm: $4 >> ${cmd_log}
echo Abundance Analysis Type: $5 >> ${cmd_log}

if [ $DEL_META = true ]
then
    rm pipeline/chira/${DATA}*
    rm pipeline/find_deletion/ALL_output/${DATA}*
    rm pipeline/predict_site/scan_output/${DATA}*
    rm pipeline/predict_site/mir_output/${DATA}*
    rm pipeline/predict_site/up_output/${DATA}*
    rm pipeline/data_processing/after_preprocess/${DATA}*
    rm pipeline/${OUTPUT}
fi
echo "Program complete successfully."
