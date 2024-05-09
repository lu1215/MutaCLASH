cp $1 hyb_file.fastq
cp $2 reg_file.fasta
cp $3 tran_file.fasta

# preprocess
. ../../preprocess.conf
suite=../../bio_tool/clash_analyst
make -f ${suite}/bin/makefile preprocess qc=trim_galore trim=${trim} link=${link} len=${len} slen=${slen} rc=${rc} fd=${fd} in=hyb_file.fastq

# build reference and run pipeline
if [ $4 = "pir" ]
then
    make -f ${suite}/bin/makefile build reg=reg_file.fasta tran=tran_file.fasta in=hyb_file.fastq
    make -f ${suite}/bin/makefile detect way=pir llen=17 reg_mis=0 tran_mis=0 hmax=10 reg=reg_file.fasta tran=tran_file.fasta in=hyb_file.fastq
    # make -f ${suite}/bin/makefile analyse reg=reg_file.fasta tran=tran_file.fasta in=hyb_file.fastq way=pir

elif [ $4 = "hyb" ]
then
    make -f ${suite}/bin/makefile build reg=reg_file.fasta tran=tran_file.fasta in=hyb_file.fastq
    make -f ${suite}/bin/makefile detect way=hyb hval=0.1 hmax=10 gmax=4 reg=reg_file.fasta tran=tran_file.fasta in=hyb_file.fastq
    # make -f ${suite}/bin/makefile analyse reg=reg_file.fasta tran=tran_file.fasta in=hyb_file.fastq way=hyb

    cut -d , -f 1,2,13,15-17,19,20 --complement hyb_file_step4.csv > output/${5}_${4}.csv
    sed -i 's/\r//g' output/${5}_${4}.csv

elif [ $4 = "clan" ]
then
    make -f ${suite}/bin/makefile build reg=reg_file.fasta tran=tran_file.fasta in=hyb_file.fastq
    make -f ${suite}/bin/makefile detect way=clan llen=17 hmax=10 gmax=4 reg=reg_file.fasta tran=tran_file.fasta in=hyb_file.fastq
    # make -f ${suite}/bin/makefile analyse reg=reg_file.fasta tran=tran_file.fasta in=hyb_file.fastq way=clan

    cut -d , -f 1,2,13,15-17,19,20 --complement hyb_file_step4.csv > output/${5}_${4}.csv
    sed -i 's/\r//g' output/${5}_${4}.csv
fi

# prepare data
python csv_to_fasta.py --type input --input hyb_file_step2.csv --output output/${5}.fa
python fasta_to_csv.py --type regulator --input $2 --output ${2%.*}.csv
python fasta_to_csv.py --type transcript --input $3 --output ${3%.*}.csv
mv hyb_file.fastq_trimming_report.txt output/${5}_trimming.log
rm hyb_file*
rm *.fasta
rm *.csv
rm -r bowtieFile
rm -r idFile
