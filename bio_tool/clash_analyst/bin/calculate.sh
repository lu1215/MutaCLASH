
result=result.tsv
printf "%-15s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s\n" 'file' "clan_gen" "clan_pair" "clan_ori" "clan_cent" "hyb_gene" "hyb_pair" "hyb_ori" "hyb_cent" "pir_gene" "pir_pair" "pir_ori" "pir_cent"> $result
for i in {87..92};do
    f=GSM12194${i}
    cd $f
    way=pir
    declare -i ${way}_gene=$(wc -l $way"_gene"|cut -d" " -f 1)
    declare -i ${way}_pair=$(wc -l $way"_pair"|cut -d" " -f 1)
    declare -i ${way}_ori=$(wc -l $way"_ori"|cut -d" " -f 1)
    pir_cent=$(echo "scale=2;${pir_pair}*100/${pir_ori}"|bc)

    way=clan
    declare -i ${way}_gene=$(wc -l $way"_gene"|cut -d" " -f 1)
    declare -i ${way}_pair=$(wc -l $way"_pair"|cut -d" " -f 1)
    declare -i ${way}_ori=$(wc -l $way"_ori"|cut -d" " -f 1)
    clan_cent=$(echo "scale=2;${clan_pair} * 100 / ${clan_ori}"|bc)
    way=hyb
    declare -i ${way}_gene=$(wc -l $way"_gene"|cut -d" " -f 1)
    declare -i ${way}_pair=$(wc -l $way"_pair"|cut -d" " -f 1)
    declare -i ${way}_ori=$(wc -l $way"_ori"|cut -d" " -f 1)
    hyb_cent=$(echo "scale=2;${hyb_pair} * 100 / ${hyb_ori}"|bc)
    cd ..
    printf "%-15s %10s %10s %10s %10.2f %10s %10s %10s %10.2f %10s %10s %10s %10.2f\n" $f $clan_gene $clan_pair $clan_ori $clan_cent $hyb_gene $hyb_pair $hyb_ori $hyb_cent $pir_gene $pir_pair $pir_ori $pir_cent>> $result

done



way=clan
declare -i ${way}_gene=$(cat */$way"_gene"|sort -u|wc -l|cut -d" " -f 1)
declare -i ${way}_pair=$(cat */$way"_pair"|wc -l|cut -d" " -f 1)
declare -i ${way}_ori=$(cat */$way"_ori"|wc -l|cut -d" " -f 1)
clan_cent=$(echo "scale=2;${clan_pair} * 100 / ${clan_ori}"|bc)

way=hyb
declare -i ${way}_gene=$(cat */$way"_gene"|sort -u|wc -l|cut -d" " -f 1)
declare -i ${way}_pair=$(cat */$way"_pair"|wc -l|cut -d" " -f 1)
declare -i ${way}_ori=$(cat */$way"_ori"|wc -l|cut -d" " -f 1)
hyb_cent=$(echo "scale=2;${hyb_pair} * 100 / ${hyb_ori}"|bc)

way=pir
declare -i ${way}_gene=$(cat */$way"_gene"|sort -u|wc -l|cut -d" " -f 1)
declare -i ${way}_pair=$(cat */$way"_pair"|wc -l|cut -d" " -f 1)
declare -i ${way}_ori=$(cat */$way"_ori"|wc -l|cut -d" " -f 1)
pir_cent=$(echo "scale=2;${pir_pair}*100/${pir_ori}"|bc)


printf "%-15s %10s %10s %10s %10.2f %10s %10s %10s %10.2f %10s %10s %10s %10.2f\n" "ALL" $clan_gene $clan_pair $clan_ori $clan_cent $hyb_gene $hyb_pair $hyb_ori $hyb_cent $pir_gene $pir_pair $pir_ori $pir_cent>> $result
