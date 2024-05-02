function FindPair(){
    f1=$1
    way=$2
    miRNA=$3
    mRNA=$4
    cut -d$'\t' -f2,3 $f1/$way".tab"|grep "|"$mRNA"|" |grep $miRNA >>  "real_pair/"$way"_"$f1
    find=$?
    if [ $find -eq 0 ];then
        echo $miRNA,$mRNA >> "real_pair/"$way"_"$f1"_result"
    fi

}

for i in {87..87};do
    f="GSM12194"$i
    echo $f
    printf "" >  "real_pair/clan_"$f
    printf "" >  "real_pair/clan_"$f"_result"
    printf "" >  "real_pair/hyb_"$f
    printf "" >  "real_pair/hyb_"$f"_result"
    printf "" >  "real_pair/pir_"$f
    printf "" >  "real_pair/pir_"$f"_result"
    while read line
    do
        miRNA=$(cut -d, -f1 <<< $line)
        mRNA=$(cut -d, -f2 <<< $line)
        FindPair $f clan $miRNA $mRNA
        FindPair $f hyb $miRNA $mRNA
        FindPair $f pir $miRNA $mRNA

    done <hsa_MTI.csv
done
