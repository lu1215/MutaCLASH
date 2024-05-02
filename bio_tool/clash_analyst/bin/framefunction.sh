

function merge_csv(){
    local file1=$1
    local file2=$2
    local file1_col=$3
    local file2_col=$4
    local outfile=$5
    local merge_head=$6

    if [ $merge_head -eq 1 ];then
        local header=$(paste -d, \
            <(head -n 1 $file1) \
            <(head -n 1 $file2|cut -d, -f $file2_col --complement))

        LC_ALL=C join -t, -1 $file1_col -2 $file2_col \
            <(sed '1d' $file1|LC_ALL=C sort -t, -k $file1_col,$file1_col)\
            <(sed '1d' $file2|LC_ALL=C sort -t, -k $file2_col,$file2_col)\
            | sed "1i $header"\
            > $outfile
    else
        LC_ALL=C join -t, -1 $file1_col -2 $file2_col \
            <(LC_ALL=C sort -t, -k $file1_col,$file1_col $file1)\
            <(LC_ALL=C sort -t, -k $file2_col,$file2_col $file2)\
            > $outfile
    fi
}



function find_col(){
    #echo "find_col:$1 from $2"

    local search_col=$1
    local filename=$2

    local col=$(head -n 1 $filename|sed 's/ /_/g'|awk -F, -v search_col=$search_col 'BEGIN{col=0}{for(i=1;i<=NF;i++){if($i==search_col){col=i}};print col}')
    echo $col

    #if [ $col -eq 0 ]
    #then
        #echo "You don't have \"$search_col\" column in $filename"
        #exit 1
    #fi
    #return $col
}


# call arguments verbatim:
$@
