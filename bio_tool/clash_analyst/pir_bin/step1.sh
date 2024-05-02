#!/bin/bash
#set -u -e
# example
#time bash step1.sh -i SRR6512653.fastq -o pra1.csv -l 17 -a AGATCGGAAGAG

shell_folder=$(cd "$(dirname "$0")";pwd)
project_folder=$(dirname $shell_folder)
. ${shell_folder}/framefunction.sh



# check whether use options or not
# ---------------------------------------
if [ $# -lt 1 ];then
	echo "Try '$0 -h' for more information."
	exit 1
fi

# default parameter
#-----------------------------------
length=17
output="step1.csv"
phred_score=20
use_trim=1


# -h option
#-----------------------------------
function usage(){
cat << -EOF-
Usage: 
    $0 -i <input> -o <outfile> -l <length> -a <adaptor> -q <Phred Qual> 
Options:
    -h display this help and exit

    -i input file(fastq)

    -o output file(csv)
       (sequence,read count)
       default step1.csv

    -l length(greater equal)
       default 17

    -a Adaptor for trim_galore
       If not specified,trim_galore will auto detect.

   -q Phred Qual(default 20)
-EOF-
exit 1
}


# get options
#-----------------------------------
while getopts ":i:o:l:a:q:u:h" opt
do
	case $opt in
		h)
			usage
			;;
        i)
            input="$(cd $(dirname $OPTARG);pwd)/$(basename $OPTARG)"
            ;;
        o)
            output="$(cd $(dirname $OPTARG);pwd)/$(basename $OPTARG)"
            ;;
        l)
            length=$OPTARG
            ;;
        a)
            trimmed_seq=$OPTARG
            ;;
        q)
            phred_score=$OPTARG
            ;;
        u)
            use_trim=$OPTARG
            ;;

		*)
			echo -e "$0: invalid option -- 'x'\nTry '$0 -h' for more information."
            exit 1
			;;
	esac
done

#if [ ! $trimmed_seq ];then
    #echo "you need to input '-a <adaptor>'"
    #exit 1
#fi

if [ ! $input ];then
    echo "you need to input '-i <input>'"
    exit 1
fi

temp_path=$(dirname $output)

# check file/directory exist
checkFile $input

# check the \n format,if file from dos
# use dos2unix
checkNewline $input



echo ----------parameter-----------
echo input=$input
echo output=$output
echo length=$length
echo trimmed_seq=$trimmed_seq
echo ------------------------------


# for trim_galore output file name
trimmed_input=$(basename ${input%.*})"_trimmed.fq"
clipped_input=$(basename ${input%.*})"_clipped_qf.fastq"
clipped_tab=$(basename ${input%.*})"_clipped_qf.tab"
comp_fasta=$(basename ${input%.*})"_comp.fasta"
echo $trimmed_input



if test "$use_trim" == "1";then

    if test "$trimmed_seq" == "" -o "$trimmed_seq" == "None";then

        trim_galore --length $length --dont_gzip -o $temp_path -q $phred_score $input 

    else

        trim_galore --length $length --dont_gzip -a $trimmed_seq -o $temp_path -q $phred_score $input

    fi
else
    echo "Not use trim_galore"
    cp $input ${temp_path}"/"$trimmed_input
fi

awk -f ${shell_folder}/sequence.awk ${temp_path}"/"$trimmed_input |sort| uniq -c |awk -v OFS="," 'BEGIN{print "sequence,read_count"}{print $2,$1}' > $output

# for hyb
mv ${temp_path}"/"$trimmed_input ${temp_path}"/"$clipped_input 
touch ${temp_path}"/"${clipped_tab}
sed '1d' $output|awk -F, '{print ">"NR"_"$2;print $1}' > ${temp_path}"/"${comp_fasta} 


