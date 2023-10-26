# check input file
INPUT_NAME=${1%.*}
INPUT_EX=${1##*.}

if [ $INPUT_EX = "csv" ]; then
    python csv_to_fasta.py --input $1 --output ${INPUT_NAME}.fa --type input
elif [ $INPUT_EX = "fasta" ]; then
    cp $1 ${INPUT_NAME}.fa
fi

# check regulator file
REG_NAME=${2%.*}
REG_EX=${2##*.}

if [ $REG_EX = "csv" ]; then
    python csv_to_fasta.py --input $2 --output ${REG_NAME}.fa --type regulator
elif [ $REG_EX = "fasta" ]; then
    python fasta_to_csv.py --input $2 --output ${REG_NAME}.csv --type regulator
    cp $2 ${REG_NAME}.fa
elif [ $REG_EX = "fa" ]; then
    python fasta_to_csv.py --input $2 --output ${REG_NAME}.csv --type regulator
fi

# check transcript file
TRANS_NAME=${3%.*}
TRANS_EX=${3##*.}
if [ $TRANS_EX = "csv" ]; then
    python csv_to_fasta.py --input $3 --output ${TRANS_NAME}.fa --type transcript
elif [ $TRANS_EX = "fasta" ]; then
    python fasta_to_csv.py --input $3 --output ${TRANS_NAME}.csv --type transcript
    cp $3 ${TRANS_NAME}.fa
elif [ $TRANS_EX = "fa" ]; then
    python fasta_to_csv.py --input $3 --output ${TRANS_NAME}.csv --type transcript
fi
