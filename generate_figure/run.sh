TOOL=$3
TYPE=$4
ARGS="--basename $1 --inputname $2 --tool $3 --norm_factor $5"

python rm_fig.py

echo "basic info"
python plot_nofilter_basic.py ${ARGS} > log/${1}_info.log

echo "plot CIMS"
if [ $TOOL = "pirScan" ]; then
    python plot_nofilter_pirScan.py ${ARGS} > log/${1}_pirScan.log
elif [ $TOOL = "miRnada" ]; then
    python plot_nofilter_miRanda.py ${ARGS} > log/${1}_miRnada.log
elif [ $TOOL = "RNAup" ]; then
    python plot_nofilter_RNAup.py ${ARGS} > log/${1}_RNAup.log
fi

echo "plot enrichment"
if [ $TYPE = "abu" ]; then
    python plot_nofilter_mRNA_abu.py ${ARGS} > log/${1}_enrichment.log
else
    python plot_nofilter_22G_abu.py ${ARGS} > log/${1}_enrichment.log
fi
