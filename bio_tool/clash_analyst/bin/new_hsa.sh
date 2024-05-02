for i in {87..92};do
    f="GSM12194"$i
    mkdir real_pair/$f
    echo $f
    way=pir
    cat $f/${way}.tab|awk -F$'\t' 'BEGIN{OFS="\t"}{split($2,miRNA," ");split($3,mRNA,"|");print miRNA[1],mRNA[6]}'> real_pair/$f/${way}"_ori"
    comm -12 u_hsa_MTI.tab <(sort -u real_pair/$f/${way}"_ori")> real_pair/$f/$way"_gene"
    cat real_pair/$f/$way|xargs -I{} grep "^{}$" real_pair/$f/$way"_ori" > real_pair/$f/$way"_pair"
    way=hyb
    cat $f/${way}.tab|awk -F$'\t' 'BEGIN{OFS="\t"}{split($2,miRNA," ");split($3,mRNA,"|");print miRNA[1],mRNA[6]}'> real_pair/$f/${way}"_ori"
    comm -12 u_hsa_MTI.tab <(sort -u real_pair/$f/${way}"_ori")> real_pair/$f/$way"_gene"
    cat real_pair/$f/$way|xargs -I{} grep "^{}$" real_pair/$f/$way"_ori" > real_pair/$f/$way"_pair"
    way=clan
    cat $f/${way}.tab|awk -F$'\t' 'BEGIN{OFS="\t"}{split($2,miRNA," ");split($3,mRNA,"|");print miRNA[1],mRNA[6]}'> real_pair/$f/${way}"_ori"
    comm -12 u_hsa_MTI.tab <(sort -u real_pair/$f/${way}"_ori")> real_pair/$f/$way"_gene"
    cat real_pair/$f/$way|xargs -I{} grep "^{}$" real_pair/$f/$way"_ori" > real_pair/$f/$way"_pair"
done


#f=GSM1219487
#way=clan



