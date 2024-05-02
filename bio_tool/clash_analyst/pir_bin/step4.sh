#!/bin/bash

#time bash step4.sh -i step3.csv -o step4.csv -r ~/pra/input/mRNA_sequence.csv -m 0 -b 1

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
bflag=1
outfile=step4.csv


# -h option
function usage(){
cat << -EOF-
Usage: 
    $0 -i <input> -o <output> -r <reference> -m <mismatch> -b <build>
Options:
    -h display this help and exit

    -i input file(csv)
       (need "remain0","remain_seq" column name)

    -o output file(csv)
       default step4.csv

    -r reference file(csv)
       (need "sequence" column name)

    -m mismatch count(0,1,2)
       default 0

    -b 1 use bowtie build
       0 not use bowtie build
       default 1
-EOF-
exit 1
}

# get options
# ---------------------------------------
while getopts ":i:o:r:m:s:b:h" opt
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
base_ref=$(basename $ref)
base_inp=$(basename $inp)
bowtie_path=${temp_path}/bowtieFile
id_file_in="${temp_path}/idFile/${base_inp}.csv"
id_file_ref="${temp_path}/idFile/${base_ref}.csv"
bowtie_extract=$bowtie_path/Reads_col2.bwt
merge_RNA=${temp_path}"/merge_"${base_ref}".csv"



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
csvTofasta $refile transcript sequence 


# for CPU
CPU_num
cpu_num=$?
echo -------------USE CPU number : ${cpu_num} ----------------


find_col remain0 $infile
remain_id_col=$?
find_col remain_seq $infile
remain_seq_col=$?

awk -F, -v id=$remain_id_col -v seq=$remain_seq_col 'NR==1 {next}{printf ">%s\n%s\n",$id,$seq}' $infile | sed "s/U/T/g" > ${temp_path}"/"${base_inp}"temp.fasta" 



# create a id,seq,seq_len file 
# replace sequence U to T
# at idFile/
addID $refile transcript $id_file_ref


# bowtie
echo "----------bowtie-----------------"
if [ "$bflag" = "1" ]
then
    bowtie-build --threads $cpu_num $ref"temp.fasta" "${bowtie_path}/${base_ref}.fa" > /dev/null 2>&1
else
    echo "you don't use bowtie-build"
fi

bowtie --threads $cpu_num -f -a -v $mismatch --norc ${bowtie_path}/${base_ref}".fa" ${temp_path}"/"${base_inp}"temp.fasta" $bowtie_path/Reads2.bwt > /dev/null 2>&1


# check the bowtie output is empty or not
declare -i line=0
line=$(cat $bowtie_path/Reads2.bwt|wc -l)
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
echo "remain0,transcript0,rem_tran_target_pos,rem_tran_mismatch_count" > $bowtie_extract
if [ $mismatch -gt 0 ];then
    awk -f ${shell_folder}/mismatch.awk $bowtie_path/Reads2.bwt >> $bowtie_extract
else
    awk '{printf "%s,%s,%s-%s,0\n",$1,$3,$4+1,$4+length($6)}' $bowtie_path/Reads2.bwt >> $bowtie_extract
fi


echo "---------merge file-------------"
merge_csv $id_file_ref 1 $bowtie_extract 2 $merge_RNA 


find_col remain0 $merge_RNA
RNA_remain_col=$?
echo aaa:$RNA_remain_col

merge_csv $infile $remain_id_col $merge_RNA $RNA_remain_col $outfile


rm ${temp_path}"/"${base_inp}"temp.fasta" ${ref}"temp.fasta" $merge_RNA

