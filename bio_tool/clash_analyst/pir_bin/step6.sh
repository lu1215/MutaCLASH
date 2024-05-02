#!/bin/bash

#time bash step6.sh -i step5.csv -o step6.csv -t original_name/mRNA_sequence.csv -r original_name/piRNA1000.csv -p 0 -u -10 -f -10

shell_folder=$(cd "$(dirname "$0")";pwd)
. ${shell_folder}/framefunction.sh
. ${shell_folder}/filefunction.sh



# check whether use options or not
# ---------------------------------------
if [ $# -lt 1 ];then
	echo "Try '$0 -h' for more information."
	exit 1
fi

# default parameter
#-----------------------------------
piScore=None
upScore=None
foldScore=None
outfile=step6.csv


# -h option
function usage(){
cat << -EOF-
Usage: 
    $0 -i <input> -o <output> -r <original regulator file> -t <original transcript file> -p <GU_targeting_score> -u <RNAup_score> -f <RNAfold_MFE>
Options:
    -h display this help and exit

    -i input file(csv)
       (need following column name)
       "transcript_name",
       "regulator_name",
       "hybrid_seq",
       "RNAup_score","RNAfold_MFE",

    -o output file (csv)
       same as <input> format.
       default step6.csv
       Another table output will be 
           OUTPUT_transcript_regulator.csv,
           OUTPUT_regulator_transcript.csv,
           OUTPUT_hybrid_transcript.csv

    -t original transcript_name file
       (need "transcript_name" column name)

    -r original regulator_name file
       (need "regulator_name" column name)

    -p select "GU_targeting_score" (greater equal)
       (if you want use this option,you need to have "GU_targeting_score" column name)
       default None

    -u select "RNAup_score" (less equal)
       default None

    -f select "RNAfold_MFE" (less equal)
       default None


-EOF-
exit 1
}


# ---------------------------------------
while getopts ":i:o:t:r:p:f:u:h" opt
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
        t)
            tranfile="$(cd $(dirname $OPTARG);pwd)/$(basename $OPTARG)"
            ;;
        r)
            regfile="$(cd $(dirname $OPTARG);pwd)/$(basename $OPTARG)"
            ;;
        p)
            piScore=$OPTARG
            ;;
        
        f)
            foldScore=$OPTARG
            ;;
        u)
            upScore=$OPTARG
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
if [ ! $regfile ];then
    echo "you need to input '-r <original regulator file>'"
    exit 1
fi

if [ ! $tranfile ];then
    echo "you need to input '-t <original transcript file>'"
    exit 1
fi

# check file/directory exist
checkFile $infile
checkFile $regfile
checkFile $tranfile

# check the \n format,if file from dos
# use dos2unix
checkNewline $infile
checkNewline $regfile
checkNewline $tranfile

temp_path=$(dirname $outfile)


if [ "$piScore" != "None" ];then
    echo "choose piScore >= $piScore"
    find_col GU_targeting_score $infile
    p_col=$?
else 
    piScore=-1000
    p_col=-1
fi


if [ "$foldScore" != "None" ];then
    echo "choose RNAfoldScore <= $foldScore"
else
    foldScore=1000
fi

if [ "$upScore" != "None" ];then
    echo "choose RNAupScore <= $upScore"
else 
    upScore=1000
fi




find_col transcript_name $infile
RNA_col=$?
find_col regulator_name $infile
input_col=$?
find_col hybrid_seq $infile
reference_col=$?
find_col RNAfold_MFE $infile
fold_col=$?
find_col RNAup_score $infile
up_col=$?

find_col transcript_name $tranfile
tran_col=$?
find_col regulator_name $regfile
reg_col=$?


base_inp=$(basename ${infile%.*})
out=${outfile%.*}

if [ $p_col -eq -1 ];then
    awk -F, -v upScore=$upScore -v foldScore=$foldScore -v fold_col=$fold_col -v up_col=$up_col '($up_col <= upScore && $fold_col <= foldScore)||NR==1 {$1=$1;print $0}' OFS="," $infile > $outfile
else
    awk -F, -v piScore=$piScore -v upScore=$upScore -v foldScore=$foldScore -v fold_col=$fold_col -v up_col=$up_col -v p_col=$p_col '($p_col >= piScore && $up_col <= upScore && $fold_col <= foldScore)||NR==1 {$1=$1;print $0}' OFS="," $infile > $outfile
fi



cut -d, -f $input_col,$RNA_col $outfile|sed '1d'|LC_ALL=C sort -u |sed '1i new line'> ${temp_path}/t.csv
if [ $input_col -lt $RNA_col ];then
    awk -F, -v key=2 -v value=1 -f ${shell_folder}/key_value.awk ${temp_path}/t.csv|LC_ALL=C sort -k 1,1>${temp_path}/temp_transcript_regulator.csv
    awk -F, -v key=1 -v value=2 -f ${shell_folder}/key_value.awk ${temp_path}/t.csv|LC_ALL=C sort -k 1,1>${temp_path}/temp_regulator_transcript.csv
    rm ${temp_path}/t.csv
else 
    awk -F, -v key=1 -v value=2 -f ${shell_folder}/key_value.awk ${temp_path}/t.csv|LC_ALL=C sort -k 1,1>${temp_path}/temp_transcript_regulator.csv
    awk -F, -v key=2 -v value=1 -f ${shell_folder}/key_value.awk ${temp_path}/t.csv|LC_ALL=C sort -k 1,1>${temp_path}/temp_regulator_transcript.csv
    rm ${temp_path}/t.csv
fi

# hybrid to target pair
echo "hybrid_seq,target_pair_sum,target_pair_name" > ${out}_hybrid_transcript.csv
awk -F, -v hyb=$reference_col -v reg=$input_col -v tran=$RNA_col -f ${shell_folder}/hybrid_target.awk $outfile| LC_ALL=C sort -u|sed '1i new line'|awk -F, -v key=1 -v value=2 -f ${shell_folder}/key_value.awk |LC_ALL=C sort >> ${out}_hybrid_transcript.csv


echo "transcript_name,regulator_sum,regulator_name" > ${out}_transcript_regulator.csv
# join also need LC_ALL
cut -d, -f $tran_col $tranfile|sed '1d'|LC_ALL=C sort -t, -k 1,1 |LC_ALL=C join -t, -a 1 -e 0 -o 1.1,2.2,2.3 - ${temp_path}/temp_transcript_regulator.csv >> ${out}_transcript_regulator.csv

echo "regulator_name,transcript_sum,transcript_name" > ${out}_regulator_transcript.csv
cut -d, -f $reg_col $regfile|sed '1d'|LC_ALL=C sort -t, -k 1,1 |LC_ALL=C join -t, -a 1 -e 0 -o 1.1,2.2,2.3 - ${temp_path}/temp_regulator_transcript.csv >> ${out}_regulator_transcript.csv

rm ${temp_path}/temp_transcript_regulator.csv ${temp_path}/temp_regulator_transcript.csv

