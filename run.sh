# read path
READ=../../$1
# regulator path
REG=../../$2
# target path
TAR=../../$3
# data base name
DATA=$(basename ${READ})
DATA=${DATA%.*}

# binding-site & enrichment analysis
# [pirScan/miRnada/RNAup]
TOOL=$4
# [site/region/up/abu]
TYPE=$5

# remove metadatas
DEL_META=false

# ===========================

echo "Step0. chira"
cd pipeline/chira
REG=${REG%.*}.fa
TAR=${TAR%.*}.fa
# run.sh [fold_name] [regulator] [target] [thread] [seed_length(12)] [gap_penalty(6)] [mismatch_penalty(4)] [score_cutoff(18)]
sh run.sh ${DATA} ${READ} ${REG} ${TAR} 8 12 6 4 18
cd ..
OUTPUT=chira/${DATA}_extract_dir/${DATA}.csv

# --------------------------

echo "Step1. find deletion"
cd find_deletion
REG=${REG%.*}.csv
TAR=${TAR%.*}.csv
sh run.sh ../${OUTPUT} ${TAR}
cd ..
OUTPUT=find_deletion/ALL_output/${DATA}_step1.csv

# --------------------------

echo "Step2. predict site"
cd predict_site
# [n/extend_length]
EXTEND=n
sh run_pir.sh ../${OUTPUT} ${REG} ${TAR} ${EXTEND}
OUTPUT=predict_site/scan_output/${DATA}_step1_scan.csv
sh run_mir.sh ../${OUTPUT} ${REG} ${TAR} ${EXTEND}
OUTPUT=predict_site/mir_output/${DATA}_step1_scan_mir.csv
cd ..

# --------------------------

echo "Step3. RNAup"
cd RNAup
sh run.sh ../${OUTPUT} ${REG} ${TAR} ${EXTEND}
OUTPUT=RNAup/RNAup_output/${DATA}_step1_scan_mir_RNAup.csv
cd ..

# --------------------------

echo "Step4. data processing"
cd data_processing
sh run.sh ../${OUTPUT} ${TAR}
OUTPUT=data_processing/after_preprocess/${DATA}_step1_scan_mir_RNAup_final.csv
cd ..

# --------------------------

echo "Step5. enrichment"
cd induce_22g
# [n/extend_length]
EXTEND=25
sh run.sh ../${OUTPUT} ${REG} ${TAR} ${EXTEND} ${TYPE}
if [ $TYPE = "abu" ]
then
    OUTPUT=induce_22g/add_abu_info/abu_${EXTEND}_${DATA}_step1_scan_mir_RNAup_final.csv
else
    OUTPUT=induce_22g/add_22g_info/22g_${TYPE}_${EXTEND}_${DATA}_step1_scan_mir_RNAup_final.csv
fi
cd ..

# --------------------------

echo "Step6. generate figure"
cd generate_figure
# normalization factor
FACTOR=811.03
sh run.sh ${DATA} ../${OUTPUT} ${TOOL} ${TYPE} ${FACTOR}
cd ../../

# --------------------------

echo "Step7. collect files"
cp pipeline/${OUTPUT} data/output/
cp -r pipeline/generate_figure/figure data/output/

if [ $DEL_META = true ]
then
    rm pipeline/chira/${DATA}*
    rm pipeline/find_deletion/ALL_output/${DATA}*
    rm pipeline/predict_site/scan_output/${DATA}*
    rm pipeline/predict_site/mir_output/${DATA}*
    rm pipeline/RNAup/RNAup_output/${DATA}*
    rm pipeline/data_processing/after_preprocess/${DATA}*
    rm pipeline/${OUTPUT}
fi
echo "Program complete successful."
