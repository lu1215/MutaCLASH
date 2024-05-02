#!/bin/bash


function CPU_num(){
    declare -i cpu_num=$(cat /proc/cpuinfo | grep "processor" | wc -l)   
    if [ $cpu_num -ge 4 ];then
        cpu_num=$cpu_num-2
    fi
    return $cpu_num
}

function createDir(){
    if [ ! -d $1  ]; then
          mkdir -p $1
    fi
}


# check file/directory exist
# ---------------------------------------
function checkFile(){
    if [ ! -e $1 ];then
        echo "plz check your $1 whether exists or not."
        exit 1
    fi
}


# check the \n format,if file from dos
# use dos2unix
# ---------------------------------------
function checkNewline(){
    if head -n 1 $1 |cat -v - |grep -q "M$";then dos2unix $1; fi
}

# select
# ---------------------------------------
#select bigger number
function selectBigger(){
    local col_name=$1
   local count=$2
    local filename=$3
    local outputfile=$4

    find_col $col_name $filename 
    local col_num=$?

    awk -F, -v count=$count -v col=$col_num '$col>=count||NR==1 {$1=$1;print $0}' OFS="," $filename > $outputfile
}
#selectNum 5 2 /home/bba753951/master_project/clash_data/CLASH.WT.inserts.uniq.reads lll.csv


#select smaller number
function selectSmaller(){
    local col_name=$1
    local count=$2
    local filename=$3
    local outputfile=$4

    find_col $col_name $filename 
    local col_num=$?

    awk -F, -v count=$count -v col=$col_num '$col<=count||NR==1 {$1=$1;print $0}' OFS="," $filename > $outputfile
}


# ---------------------------------------
# input : search_col, filename
# return : int(index start from 1) 
function find_col(){
    echo "find_col:$1 from $2"

    local search_col=$1
    local filename=$2

    local col=$(head -n 1 $filename|sed 's/ /_/g'|awk -F, -v search_col=$search_col 'BEGIN{col=0}{for(i=1;i<=NF;i++){if($i==search_col){col=i}};print col}')

    if [ $col -eq 0 ]
    then
        echo "You don't have \"$search_col\" column in $filename"
        exit 1
    fi
    return $col
}


#find_col input0 ../output/step1/merge_inputID.csv
#echo $?


# ---------------------------------------
function find_col2(){
    echo "find_col:$1 from $2"

    local search_col=$1
    local filename=$2

    local col=$(head -n 1 $filename|sed 's/ /_/g'|awk -F, -v search_col=$search_col 'BEGIN{col=0}{for(i=1;i<=NF;i++){if($i==search_col){col=i}};print col}')
    return $col
}

# ---------------------------------------
#input:$1 (fasta)
#hyread_file:$2 (csv)
#rnafold_result:$3 (csv)
#return:$3
#note:function name can't use RNAfold!
function rnafold(){

    RNAfold -j0 -d2 --noLP --noPS --noconv < $1|
    # extract score
    grep -P -o "(?<=\s\().*(?=\)$)"|
    # remove blank
    sed -r 's/\s+//g'|
    # add column name
    sed '1i RNAfold_MFE'|
    # merge column
    paste -d, $2 - > $3 
}
#time rnafold readcount5.fasta readcount5.csv bbb

# merge two csv file
# sort before join
function merge_csv(){
    local file1=$1
    local file1_col=$2
    local file2=$3
    local file2_col=$4
    local outfile=$5

    LC_ALL=C sort -t, -k $file1_col,$file1_col -o $file1 $file1
    LC_ALL=C sort -t, -k $file2_col,$file2_col -o $file2 $file2

    LC_ALL=C join -t, -1 $file1_col -2 $file2_col $file1 $file2 > $outfile

}

