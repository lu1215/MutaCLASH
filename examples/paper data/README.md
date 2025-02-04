PRG-1 CLASH data processing
```
git clone https://github.com/lu1215/MutaCLASH.git
cd MutaCLASH/data/input
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR6512652/SRR6512652
../../bio_tool/sratoolkit/bin/fastq-dump SRR6512652
cd ../..
pip install requirements.txt
# change REGION=10/0/-15/-30(run_all.sh: 231) 
sh run_all.sh --input data/input/SRR6512652.fastq --regulator data/reference/piRNA_WS275.fa --transcript data/reference/mRNA_WS275.fa --algorithm pirScan --abundance_type site
```

ALG-1 CLASH data processing
```
git clone https://github.com/lu1215/MutaCLASH.git
cd MutaCLASH/data/input
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR3882949/SRR3882949
../../bio_tool/sratoolkit/bin/fastq-dump SRR3882949
cd ../..
pip install requirements.txt
# change REGION=200/140/100/60(run_all.sh: 231) 
sh run_all.sh --input data/input/SRR3882949.fastq --regulator data/reference/miRNA_WS275.fa --transcript data/reference/mRNA_WS275.fa --algorithm miRanda --abundance_type abu
```