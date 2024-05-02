
shell_folder=$(cd "$(dirname "$0")";pwd)


if [ ! -e "plot" ];then
    echo mkdir plot
    mkdir plot
fi

name=$1
python ${shell_folder}/compare.py select $name

echo select done
pir_out=pir.tab
hyb_out=hyb.tab
clan_out=clan.tab
hybrid_file=/home/bba753951/Django/master_project/media/uploadfile/$name/hyb_file_comp.fasta
export LC_ALL=C
sort -u $pir_out |sort -t$'\t' -k 1,1 -o $pir_out 
sort -u $hyb_out |sort -t$'\t' -k 1,1 -o $hyb_out 
sort -u $clan_out |sort -t$'\t' -k 1,1 -o $clan_out 

LC_ALL=C comm -12 $pir_out $clan_out > pir_clan.tab
LC_ALL=C comm -12 $pir_out $hyb_out > pir_hyb.tab
LC_ALL=C comm -12 $hyb_out $clan_out > hyb_clan.tab
LC_ALL=C comm -12 hyb_clan.tab $pir_out > pir_hyb_clan.tab


declare -i pir_count=$(wc -l $pir_out |cut -d" " -f 1)
declare -i hyb_count=$(wc -l $hyb_out |cut -d" " -f 1)
declare -i clan_count=$(wc -l $clan_out |cut -d" " -f 1)
declare -i pir_clan_count=$(wc -l pir_clan.tab |cut -d" " -f 1)
declare -i hyb_clan_count=$(wc -l hyb_clan.tab |cut -d" " -f 1)
declare -i pir_hyb_count=$(wc -l pir_hyb.tab |cut -d" " -f 1)
declare -i pir_hyb_clan_count=$(wc -l pir_hyb_clan.tab |cut -d" " -f 1)





declare -i hybrid_count=$(expr $(wc -l ${hybrid_file}|cut -d" " -f 1) / 2)
echo hybrid_count: $hybrid_count

declare -i pir_brid=$(cut -d$'\t' -f 1 $pir_out|sort -u|wc -l)
declare -i hyb_brid=$(cut -d$'\t' -f 1 $hyb_out|sort -u|wc -l)
declare -i clan_brid=$(cut -d$'\t' -f 1 $clan_out|sort -u|wc -l)
declare -i pir_clan_brid=$(cut -d$'\t' -f 1 pir_clan.tab|sort -u|wc -l)
declare -i pir_hyb_brid=$(cut -d$'\t' -f 1 pir_hyb.tab|sort -u|wc -l)
declare -i hyb_clan_brid=$(cut -d$'\t' -f 1 hyb_clan.tab|sort -u|wc -l)
declare -i pir_hyb_clan_brid=$(cut -d$'\t' -f 1 pir_hyb_clan.tab|sort -u|wc -l)




# percent
cent_pir=$(echo "scale=2;${pir_brid}*100/${hybrid_count}"|bc)
cent_hyb=$(echo "scale=2;${hyb_brid}*100/${hybrid_count}"|bc)
cent_clan=$(echo "scale=2;${clan_brid}*100/${hybrid_count}"|bc)
cent_pir_clan=$(echo "scale=2;${pir_clan_brid}*100/${hybrid_count}"|bc)
cent_pir_hyb=$(echo "scale=2;${pir_hyb_brid}*100/${hybrid_count}"|bc)
cent_hyb_clan=$(echo "scale=2;${hyb_clan_brid}*100/${hybrid_count}"|bc)
cent_pir_hyb_clan=$(echo "scale=2;${pir_hyb_clan_brid}*100/${hybrid_count}"|bc)

echo "File Name:"$name > result.tab
echo 'original hybrid sequences(after trimmed and unique):'$hybrid_count >> result.tab
echo "">>result.tab
printf "%-15s %15s %15s %15s\n" 'pipeline' "pair_count" "hybrid_count" "hybrid_hit(%)" >> result.tab
printf "%-15s %15s %15s %15.2f\n" 'pir' $pir_count $pir_brid $cent_pir >> result.tab
printf "%-15s %15s %15s %15.2f\n" 'hyb' $hyb_count $hyb_brid $cent_hyb>> result.tab
printf "%-15s %15s %15s %15.2f\n" 'clan' $clan_count $clan_brid $cent_clan>> result.tab
printf "%-15s %15s %15s %15.2f\n" 'pir_clan' $pir_clan_count $pir_clan_brid $cent_pir_clan>> result.tab
printf "%-15s %15s %15s %15.2f\n" 'pir_hyb' $pir_hyb_count $pir_hyb_brid $cent_pir_hyb>> result.tab
printf "%-15s %15s %15s %15.2f\n" 'hyb_clan' $hyb_clan_count $hyb_clan_brid $cent_hyb_clan>> result.tab
printf "%-15s %15s %15s %15.2f\n" 'pir_hyb_clan' $pir_hyb_clan_count $pir_hyb_clan_brid $cent_pir_hyb_clan>> result.tab

echo -e "\n\n\n" >> result.tab


echo comm done

python ${shell_folder}/compare.py draw $name 
echo draw done
python ${shell_folder}/compare.py boxplot $name> mean.txt
echo boxplot done

python ${shell_folder}/compare.py pvalue $name
echo pvalue done

