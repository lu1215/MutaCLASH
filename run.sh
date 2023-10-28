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
# run.sh [fold_name] [regulator] [target] [thread(8)] [seed_length(12)] [gap_penalty(6)] [mismatch_penalty(4)] [score_cutoff(18)]
# >>>
sh run.sh ${DATA} ${READ} ${REG} ${TAR} 8 12 6 4 18
# >>>
cd ..
OUTPUT=chira/${DATA}_extract_dir/${DATA}.csv

# --------------------------

echo "Step3. find deletion"
cd find_deletion
REG=${REG%.*}.csv
TAR=${TAR%.*}.csv
# >>>
sh run.sh ../${OUTPUT} ${TAR}
# >>>
cd ..
OUTPUT=find_deletion/ALL_output/${DATA}_step1.csv

# --------------------------

echo "Step4. predict site"
cd predict_site
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

echo "Step6. enrichment"
cd induce_22g
# [n/extend_length]
EXTEND=25
# [region/site/up/abu]
TYPE=$5
# >>>
sh run.sh ../${OUTPUT} ${REG} ${TAR} ${EXTEND} ${TYPE}
# >>>
if [ $TYPE = "abu" ]
then
    OUTPUT=induce_22g/add_abu_info/abu_${EXTEND}_${DATA}_step1_scan_mir_RNAup_final.csv
else
    OUTPUT=induce_22g/add_22g_info/22g_${TYPE}_${EXTEND}_${DATA}_step1_scan_mir_RNAup_final.csv
fi
cd ..

# --------------------------

echo "Step7. generate figure"
cd generate_figure
# [pirScan/miRanda/RNAup]
TOOL=$4
# normalization factor
FACTOR=811.03
# >>>
sh run.sh ${DATA} ../${OUTPUT} ${TOOL} ${TYPE} ${FACTOR}
# >>>
cd ../../

# --------------------------

echo "Step8. collect files"
DIR=${DATA}_$(date +%Y-%m-%d_%H-%M-%S)
mkdir data/output/${DIR}
cp pipeline/${OUTPUT} data/output/${DIR}/${DATA}.csv
cp -r pipeline/generate_figure/figure data/output/${DIR}/
cp -r pipeline/generate_figure/log data/output/${DIR}/
cmd_log=data/output/${DIR}/log/command.log
touch ${cmd_log}
echo Read File: $1 >> ${cmd_log}
echo Regulator File: $2 >> ${cmd_log}
echo Transcript File: $3 >> ${cmd_log}
echp Algorithm: $4 >> ${cmd_log}
echo Enrichment Analysis Type: $5 >> ${cmd_log}

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
echo "Program complete successful."
