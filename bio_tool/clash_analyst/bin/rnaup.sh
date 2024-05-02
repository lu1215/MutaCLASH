for i in {87..92};do
    f="GSM12194"$i
    echo $f
    for i in $(seq 5 5 15);do
        file=up_$i
        echo $i
        up=$i
        way=pir
        mkdir $file/$f
        cat $f/${way}_up_${up}.tab|awk -F$'\t' 'BEGIN{OFS="\t"}{split($1,miRNA," ");split($2,mRNA,"|");print miRNA[1],mRNA[6]}'> $file/$f/${way}"_ori"
        comm -12 u_hsa_MTI.tab <(sort -u $file/$f/${way}"_ori") > $file/$f/${way}"_gene"
        cat $file/$f/${way}"_gene"|xargs -I{} grep "^{}$" $file/$f/${way}"_ori" > $file/$f/${way}"_pair"
        way=hyb
        cat $f/${way}_up_${up}.tab|awk -F$'\t' 'BEGIN{OFS="\t"}{split($1,miRNA," ");split($2,mRNA,"|");print miRNA[1],mRNA[6]}'> $file/$f/${way}"_ori"
        comm -12 u_hsa_MTI.tab <(sort -u $file/$f/${way}"_ori") > $file/$f/${way}"_gene"
        cat $file/$f/${way}"_gene"|xargs -I{} grep "^{}$" $file/$f/${way}"_ori" > $file/$f/${way}"_pair"
        way=clan
        cat $f/${way}_up_${up}.tab|awk -F$'\t' 'BEGIN{OFS="\t"}{split($1,miRNA," ");split($2,mRNA,"|");print miRNA[1],mRNA[6]}'> $file/$f/${way}"_ori"
        comm -12 u_hsa_MTI.tab <(sort -u $file/$f/${way}"_ori") > $file/$f/${way}"_gene"
        cat $file/$f/${way}"_gene"|xargs -I{} grep "^{}$" $file/$f/${way}"_ori" > $file/$f/${way}"_pair"
    done
done

