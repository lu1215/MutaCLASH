DATA=$(basename $1)
DATA=${DATA%.*}

echo D_M
python D_M_position.py --inputname $1
echo overlap
python overlap.py --inputname ${DATA}_detail.csv
echo data process
python data_preprocess.py --inputname ${DATA}_detail_with_overlap.csv --trans $2
rm tmp/${DATA}_detail_with_overlap.csv tmp/${DATA}_detail.csv
