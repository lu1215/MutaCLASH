#!/bin/bash

shell_folder=$(cd "$(dirname "$0")";pwd)
shell_parent=$(dirname $shell_folder)
. ${shell_folder}/framefunction.sh
. ${shell_folder}/filefunction.sh

# example
#bash step5.sh -i step4.csv -o step5.csv -N -20 -g 1 -p 1




# check whether use options or not
# ---------------------------------------
if [ $# -lt 1 ];then
	echo "Try '$0 -h' for more information."
	exit 1
fi

# default parameter
#-----------------------------------
sRnaScore=None
outfile=step5.csv
GU_score=None
pflag=None


# -h option
#-----------------------------------
function usage(){
cat << -EOF-
Usage: 
    $0 -i <input> -o <output> -N <RNAup_score> -g <GU target algorithm> -p <parallel>
Options:
    -h display this help and exit

    -i input file(csv)
       (need following column name)
       rem_tran_target_pos,
       remain0,transcipt0,
       remain_len,transcript_len,regulator_len,
       remain_seq,transcript_seq,regulator_seq

    -o output file(csv)
       default step5.csv

    -N select RNAup_score(less equal)
       deafult None

    -g use "GU target algorithm" or not
       only for "regulator = piRNA" or "sequence length=21"
       default None 
       if not None,you need to give a score ,then you will use "GU targeting algorithm"
       And select the 
       if you want to use this option,make sure you have already installed "python3"
       (inclue "numpy" and "pandas")

    -p parallel option
       default None
       1 for use
       if you want to use this option,make sure you have already installed "parellel"
       sudo apt-get install parallel

-EOF-
exit 1
}

# get options
# ---------------------------------------
while getopts ":i:o:N:g:p:h" opt
do
	case $opt in
		h)
			usage
			;;
        i)
            infile="$(cd $(dirname $OPTARG);pwd)/$(basename $OPTARG)"
            ;;
        o)
            outfile="$(cd $(dirname $OPTARG);pwd)/$(basename $OPTARG)"
            ;;
        N)
            sRnaScore=$OPTARG
            ;;
        g)
            GU_score=$OPTARG
            ;;
        p)
            pflag=$OPTARG
            ;;
		*)
			echo -e "$0: invalid option -- 'x'\nTry '$0 -h' for more information."
            exit 1
			;;
	esac
done

if [ ! $infile ];then
    echo "you need to input '-i <input>'"
    exit 1
fi

# check file/directory exist
checkFile $infile

# check the \n format,if file from dos
# use dos2unix
checkNewline $infile


temp_path=$(dirname $outfile)

find_col transcript_seq $infile 
RNA_seq=$?

find_col regulator_seq $infile
input_seq=$?

find_col remain_seq $infile
remain_seq=$?

find_col rem_tran_target_pos $infile
remain_RNA_pos=$?
declare -i RNApos_col=$remain_RNA_pos

find_col transcript_len $infile
RNA_len=$?

find_col regulator_len $infile
input_len=$?

find_col remain_len $infile
remain_len=$?

find_col transcript0 $infile
RNA_col=$?


find_col remain0 $infile
remain_col=$?


#echo $RNA_seq,$input_seq,$remain_seq,$remain_RNA_pos,$RNA_len,$input_len,$remain_len
echo ------- rnaup start ---------------

# add new cloumn name to output
head -n 1 $infile | sed 's/$/,RNAup_score,RNAup_pos,RNAup_target_seq,RNAup_input_seq/' > ${temp_path}/temp1.csv
# remove first row
#sed -i '1d' $infile


if [ "$pflag" = "1" ];then
    echo --------- use parallel -----------------
    cat $infile|sed '1d'|parallel --pipe --block 0.5M  awk -v RNA_seq=$RNA_seq -v input_seq=$input_seq -v remain_seq=$remain_seq -v remain_RNA_pos=$remain_RNA_pos -v RNA_len=$RNA_len -v input_len=$input_len -v remain_len=$remain_len -f ${shell_folder}/rnaup.awk  >> ${temp_path}/temp1.csv 2>/dev/null
else
    echo ---------not use parallel -----------------
    cat $infile|sed '1d'|awk -v RNA_seq=$RNA_seq -v input_seq=$input_seq -v remain_seq=$remain_seq -v remain_RNA_pos=$remain_RNA_pos -v RNA_len=$RNA_len -v input_len=$input_len -v remain_len=$remain_len -f ${shell_folder}/rnaup.awk  >> ${temp_path}/temp1.csv 2>/dev/null
fi



if [ "$sRnaScore" != "None" ];then
    echo "choose RNAup_score <= $sRnaScore"
    selectSmaller RNAup_score $sRnaScore ${temp_path}/temp1.csv ${temp_path}/temp.csv
    rm ${temp_path}/temp1.csv
else
    echo "choose RNAup_score <= $sRnaScore"
    mv ${temp_path}/temp1.csv ${temp_path}/temp.csv
fi


if [ "$GU_score" != "None" ];then

    echo ----------- use GU -------------
    time python3 ${shell_parent}/GU_targeting_algorithm/self_gu.py ${temp_path}/temp.csv 

    echo ---------- remove transcript sequence -----------------
    LC_ALL=C sort -t, -k ${remain_col}.7n -k ${RNA_col}.11n -k ${RNApos_col}n  ${temp_path}/temp_gu.csv|cut -d, -f $RNA_seq --complement > ${temp_path}/temp_gu1.csv

    echo "choose GU_targeting_score >= $GU_score"
    selectBigger GU_targeting_score $GU_score ${temp_path}/temp_gu1.csv $outfile
    rm ${temp_path}/temp.csv ${temp_path}/temp_gu.csv ${temp_path}/temp_gu1.csv

else 
    echo ----------- not use GU -------------

    LC_ALL=C sort -t, -k ${remain_col}.7n -k ${RNA_col}.11n -k ${RNApos_col}n  ${temp_path}/temp.csv|cut -d, -f $RNA_seq --complement > $outfile
    rm ${temp_path}/temp.csv
fi

