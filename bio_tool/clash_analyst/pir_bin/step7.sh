#!/bin/bash


shell_folder=$(cd "$(dirname "$0")";pwd)
. ${shell_folder}/framefunction.sh
. ${shell_folder}/filefunction.sh

#time bash step7.sh -i step6.csv -o step7.csv -G gene_transcript.csv

# check whether use options or not
# ---------------------------------------
if [ $# -lt 1 ];then
	echo "Try '$0 -h' for more information."
	exit 1
fi

# default parameter
#-----------------------------------
outfile=step7.csv


# -h option
function usage(){
cat << -EOF-
Usage: 
    $0 -i <input> -o <output> -G <Gene_transcript file>
Options:
    -h display this help and exit

    -i input file(csv)
       (need "transcript_name","regulator_name" column name)

    -o output file (csv)
       default step7.csv

    -G Gene to transcript_name file
       (need "Gene_name","transcript_name" column name)


-EOF-
exit 1
}


# ---------------------------------------
while getopts ":i:o:G:h" opt
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
        G)
            genefile="$(cd $(dirname $OPTARG);pwd)/$(basename $OPTARG)"
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
if [ ! $genefile ];then
    echo "you need to input '-G <Gene_transcript file>'"
    exit 1
fi


# check file/directory exist
checkFile $infile
checkFile $genefile

# check the \n format,if file from dos
# use dos2unix
checkNewline $infile
checkNewline $genefile


temp_path=$(dirname $outfile)

find_col Gene_name $genefile
source_geneID_col=$?

find_col transcript_name $genefile
source_RNA_col=$?


declare -i line=0
line=$(cat $infile|wc -l)
echo ------------------line $line

if [ $line -le 1 ];then
    echo ------empty file ---------
    exit 1
fi


find_col regulator_name $infile
regulator_col=$?
find_col transcript_name $infile 
transcript_col=$?

awk -F, -v key=$source_geneID_col -v value=$source_RNA_col  -f ${shell_folder}/key_value.awk $genefile|LC_ALL=C sort>${temp_path}/temp_gene_RNA.csv


if [ $regulator_col -lt $transcript_col ];then
    cut -d, -f $regulator_col,$transcript_col $infile |LC_ALL=C sort -u|awk -F, -v key=2 -v value=1 -f ${shell_folder}/key_value.awk |LC_ALL=C sort>${temp_path}/temp_RNA_regultor.csv
else 
    cut -d, -f $regulator_col,$transcript_col $infile |LC_ALL=C sort -u|awk -F, -v key=1 -v value=2 -f ${shell_folder}/key_value.awk |LC_ALL=C sort>${temp_path}/temp_RNA_regultor.csv
fi
awk -F, -f ${shell_folder}/gene_tran_reg.awk ${temp_path}/temp_RNA_regultor.csv ${temp_path}/temp_gene_RNA.csv > $outfile


rm ${temp_path}/temp_RNA_regultor.csv ${temp_path}/temp_gene_RNA.csv

