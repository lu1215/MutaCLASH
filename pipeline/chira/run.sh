# fold name
F=$1
# read
f0=$2
# regulator
f1=$3
# target
f2=$4
# single or chimeras
hybrid=$5

mkdir ${F}_map_dir
mkdir ${F}_merge_dir
mkdir ${F}_quantify_dir
mkdir ${F}_extract_dir

echo map
chira/chira_map.py -i ${f0} -o ${F}_map_dir/ -b -f1 ${f1} -f2 ${f2} -p ${6} -l1 ${7} -go1 ${8} -mm1 ${9} -s1 ${10}

echo merge
chira/chira_merge.py -b ${F}_map_dir/sorted.bed -o ${F}_merge_dir/ -f1 ${f1} -f2 ${f2}

echo quantify
chira/chira_quantify.py -b ${F}_merge_dir/segments.bed -m ${F}_merge_dir/merged.bed -o ${F}_quantify_dir/

echo extract
chira/chira_extract.py -l ${F}_quantify_dir/loci.counts -o ${F}_extract_dir/ -f1 ${f1} -f2 ${f2}

if [ $hybrid = 'single' ]; then
    echo singletons to csv
    python chira/chira_singleprocess.py -i ${f0} -t ${f2} -s ${F}_extract_dir/singletons -o ${F}_extract_dir/${F}_chira_single.csv

elif [ $hybrid = 'chimeras' ]; then
    echo chimeras to csv
    python chira/chira_chimeprocess.py -i1 ${f0} -i2 ${f1} -c ${F}_extract_dir/chimeras -o ${F}_extract_dir/${F}_chira_chimeras.csv
else
    echo unknown type of reads
fi
