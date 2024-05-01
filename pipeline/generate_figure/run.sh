TOOL=$3
TYPE=$4
ARGS="--basename $1 --inputname $2 --tool $3 --norm_factor $5 --trans $6 --figure $7"
if [ $# -ge 8 ]; then
    ARGS="${ARGS} --region $8"
fi

python rm_fig.py

echo "basic info"
python plot_nofilter_basic.py ${ARGS} > log/${1}_info.log

echo "plot CIMS"
if [ $TOOL = "pirScan" ]; then
    python plot_nofilter_pirScan.py ${ARGS} > log/${1}_pirScan.log
elif [ $TOOL = "miRanda" ]; then
    python plot_nofilter_miRanda.py ${ARGS} > log/${1}_miRanda.log
elif [ $TOOL = "RNAup" ]; then
    python plot_nofilter_RNAup.py ${ARGS} > log/${1}_RNAup.log
fi

echo "plot abundance"
if [ $TYPE = "abu" ]; then
    python plot_nofilter_mRNA_abu.py ${ARGS} > log/${1}_abundance.log
elif [ $TYPE = "region" ] || [ $TYPE = "site" ] || [ $TYPE = "up" ]; then
    python plot_nofilter_22G_abu.py ${ARGS} > log/${1}_abundance.log
fi
