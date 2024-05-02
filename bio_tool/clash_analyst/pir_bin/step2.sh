#!/bin/bash


shell_folder=$(cd "$(dirname "$0")";pwd)
. ${shell_folder}/framefunction.sh
. ${shell_folder}/filefunction.sh

# example
#time bash step2.sh -i pra1.csv -o step2.csv -T 5 -M -20



# check whether use options or not
# ---------------------------------------
if [ $# -lt 1 ];then
	echo "Try '$0 -h' for more information."
	exit 1
fi

# default parameter
#-----------------------------------
sReadcount=None
sRnaScore=None
outfile=step2.csv


# -h option
#-----------------------------------
function usage(){
cat << -EOF-
Usage: 
    $0 -i <input> -o <output> -T <read count> -M <RNAfold_MFE> 
Options:
    -h display this help and exit

    -i input file(csv)
       (need "sequence" column name)

    -o output file(csv)
       (will add a "RNAfold_MFE" column)
       default step2.csv

    -T select read count(greater equal)
       default None
       (must have "read count" or "read_count" column name)

    -M select RNAfold_MFE(less equal)
       deafult None
-EOF-
exit 1
}

# get options
# ---------------------------------------
while getopts ":i:o:T:M:h" opt
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
        T)
            sReadcount=$OPTARG
            ;;
        M)
            sRnaScore=$OPTARG
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


infile_sc=${infile%.*}"_Src.csv"
infile_fold=${infile%.*}"_fold.csv"
temp_path=$(dirname $outfile)
# select Read count
# ---------------------------------------
if [ "$sReadcount" != "None" ];then
    echo "choose read count >= $sReadcount"
    selectBigger read_count $sReadcount $infile $infile_sc
else
    echo "not choose read count"
    cp $infile $infile_sc
fi



csvTofasta $infile_sc RNAfold sequence

# RNAfold
echo "----------RNAfold-----------------"
inp=${infile_sc%.*}

rnafold $inp"temp.fasta" $infile_sc $infile_fold

# select RNAfold_MFE score
if [ "$sRnaScore" != "None" ];then
    echo "choose RNAfold_MFE <= $sRnaScore"
    selectSmaller RNAfold_MFE $sRnaScore $infile_fold $outfile
    rm $infile_fold
else
    echo "not choose RNAfold_MFE"
    mv $infile_fold $outfile
fi

rm $infile_sc $inp"temp.fasta"
echo step2 done
