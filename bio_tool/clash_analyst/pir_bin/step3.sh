#!/bin/bash

#time bash step3.sh -i ~/pra/input/piRNA1000.csv -o step3.csv -r step2.csv -m 0 -L 15 -b 1

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
mismatch=0
sRemainLen=0
bflag=1
outfile=step3.csv


# -h option
function usage(){
cat << -EOF-
Usage: 
    $0 -i <input> -o <output> -r <reference> -m <mismatch> -L <Length> -b <build>
Options:
    -h display this help and exit

    -i input/regulator file(csv)
       (need "sequence" column name)

    -o output file(csv)
       default step3.csv

    -r reference/hybrid file(csv)
       (need "sequence" column name)

    -m mismatch count(0,1,2)
       default 0

    -L select remaining sequence (greater then)
       default 0

    -b 1 use bowtie build
       0 not use bowtie build
       default 1
-EOF-
exit 1
}

# get options
# ---------------------------------------
while getopts ":i:o:r:m:L:b:h" opt
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
        r)
            refile="$(cd $(dirname $OPTARG);pwd)/$(basename $OPTARG)"
            ;;
        m)
            mismatch=$OPTARG
            ;;
        L)
            sRemainLen=$OPTARG
            ;;
		b)
            bflag=$OPTARG
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

if [ ! $refile ];then
    echo "you need to input '-r <reference>'"
    exit 1
fi

if [ $mismatch -gt 2 ];then
    echo "you just can use mismatch 0,1,2"
    exit 1
fi

# path/file
# ---------------------------------------
temp_path=$(dirname $outfile)
ref=${refile%.*}
inp=${infile%.*}
# awk can't not read file name inclue "="
base_ref=$(basename $ref|sed 's/=/_/g')
base_inp=$(basename $inp|sed 's/=/_/g')
bowtie_path=${temp_path}/bowtieFile
id_file_in="${temp_path}/idFile/${base_inp}.csv"
id_file_ref="${temp_path}/idFile/${base_ref}.csv"
bowtie_extract=$bowtie_path/Reads_col.bwt



# check file/directory exist
checkFile $infile
checkFile $refile

# check the \n format,if file from dos
# use dos2unix
checkNewline $infile
checkNewline $refile

# check dir exist or not 
# if not,create one
createDir ${temp_path}/bowtieFile
createDir ${temp_path}/idFile

# transfer csv file to fasta file
# at same Dir
csvTofasta $infile regulator sequence
csvTofasta $refile hybrid sequence 


# create a id,seq,seq_len file 
# replace sequence U to T
# at idFile/
addID $infile regulator $id_file_in
addID $refile hybrid $id_file_ref

# for CPU
CPU_num
cpu_num=$?
echo -------------USE CPU number : ${cpu_num} ----------------

# bowtie
if [ "$bflag" = "1" ]
then
    bowtie-build --threads $cpu_num $ref"temp.fasta" "$bowtie_path/$base_ref.fa" > /dev/null 2>&1
else
    echo "you don't use bowtie-build"
fi

bowtie --threads $cpu_num -f -a -v $mismatch --norc $bowtie_path/$base_ref".fa" $inp"temp.fasta" $bowtie_path/Reads.bwt > /dev/null 2>&1

# check the bowtie output is empty or not
declare -i line=0
line=$(cat $bowtie_path/Reads.bwt|wc -l)
echo ------------------line $line

if [ $line -eq 0  ];then
    echo ------not match any sequence---------
    exit 1
fi


# process Reads.bwt
# ---------------------------------------
# Reads.bwt column : input_id,ref_id,pos,mismatch_count
echo "------------extract botwie output--------------------"
echo mismatch is $mismatch
echo "regulator0,hybrid0,reg_hyb_target_pos,reg_hyb_mismatch_count" > $bowtie_extract
if [ $mismatch -gt 0 ];then
    awk -f ${shell_folder}/mismatch.awk $bowtie_path/Reads.bwt >> $bowtie_extract
else
    awk '{printf "%s,%s,%s-%s,0\n",$1,$3,$4+1,$4+length($6)}' $bowtie_path/Reads.bwt >> $bowtie_extract
fi

# merge
echo "---------merge file-------------"
merge_reg="${temp_path}/merge_reg_${base_inp}.csv"
merge_hyb="${temp_path}/merge_hyb_${base_inp}.csv"

# merge regulator file 
merge_csv $id_file_in 1 $bowtie_extract 1 $merge_reg

# find hybrid col
find_col hybrid0 $merge_reg
ref_col=$?

# merge hybrid file
merge_csv $id_file_ref 1 $merge_reg $ref_col $merge_hyb


# caculate remaining seq
echo "---------caculate remaining seq-------------"
find_col hybrid_len $merge_hyb
ref_len=$?
find_col hybrid_seq $merge_hyb
ref_seq=$?
find_col reg_hyb_target_pos $merge_hyb
target_pos=$?

awk -v ref_seq=$ref_seq -v ref_len=$ref_len -v len=$sRemainLen -v target_pos=$target_pos -f ${shell_folder}/remain.awk $merge_hyb > $outfile

rm $merge_reg $merge_hyb $inp"temp.fasta" $ref"temp.fasta"


