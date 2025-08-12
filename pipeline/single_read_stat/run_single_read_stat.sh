###################
## get parameter ##
###################

# file_name
DATA=$1

# regulator path
REG=$2

# target path
TAR=$3

# path of hybrid read analysis metadata
HYBRID_META=$4

# path of preprocessed fa


#####################################################################################
## Call sh to processing "singleton" data to get the result of single read anlysis ##
#####################################################################################
cd ../chira 

# >>
python chira/chira_singleprocess.py -i ../preprocess/output/${DATA}.fa -t ${TAR} -s ${DATA}_extract_dir/singletons -o ${DATA}_extract_dir/${DATA}_chira_single.csv
# >>

BWA_OUTPUT=chira/${DATA}_map_dir/sorted.bam
OUTPUT=chira/${DATA}_extract_dir/${DATA}_chira_single.csv

# back to single read stat folder
cd -

#####################################
## Processing mutation information ##
#####################################
cd ../find_deletion

# >>
sh run.sh chira_single ../${BWA_OUTPUT} ../${OUTPUT} ${REG} ${TAR} false
# >>

OUTPUT=find_deletion/ALL_output/${DATA}_chira_single_step1.csv

cd -

#####################################################
## calculate absolute location of mutation in mRNA ##
#####################################################

cd ../data_processing

echo D_M
# >>
python D_M_position.py --inputname ../${OUTPUT}
# >>

OUTPUT=data_processing/tmp/${DATA}_chira_single_step1_detail.csv

cd -

##########################
## single read analysis ##
##########################
# cd single_read_stat
# construct table
# >>
python construct_single_read_MUT_table.py ../${OUTPUT} output/
# >>

OUTPUT_D=output/${DATA}_chira_single_step1_detail_D_nor_rc.csv
OUTPUT_M=output/${DATA}_chira_single_step1_detail_M_nor_rc.csv

# calculate stat significance
TAR=${TAR%.*}.csv
# >>
python cal_stat.py --target ${TAR} --hybrid_csv ../${HYBRID_META} --single_D ${OUTPUT_D} --single_M ${OUTPUT_M} --output output/
# >>