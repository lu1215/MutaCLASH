# fold name
F=$1
# read
f0=$2
# regulator
f1=$3
# target
f2=$4

mkdir ${F}_map_dir
mkdir ${F}_merge_dir
mkdir ${F}_quantify_dir
mkdir ${F}_extract_dir

echo map
chira/chira_map.py -i ${f0} -o ${F}_map_dir/ -b -f1 ${f1} -f2 ${f2} -p $5 -l1 $6 -go1 $7 -mm1 $8 -s1 $9

echo merge
chira/chira_merge.py -b ${F}_map_dir/sorted.bed -o ${F}_merge_dir/ -f1 ${f1} -f2 ${f2}

echo quantify
chira/chira_quantify.py -b ${F}_merge_dir/segments.bed -m ${F}_merge_dir/merged.bed -o ${F}_quantify_dir/

echo extract
chira/chira_extract.py -l ${F}_quantify_dir/loci.counts -o ${F}_extract_dir/ -f1 ${f1} -f2 ${f2}

echo chimeras to csv
python chira/chira_postprocess.py -i1 ${f0} -i2 ${f1} -c ${F}_extract_dir/chimeras -o ${F}_extract_dir/${F}.csv
