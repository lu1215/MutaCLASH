## Instruction
- Regulator file: **piRNA_WS275.fa**
- Transcript file: **mRNA_WS275.fa**
- Input file: **C.elegans PRG-1 CLASH data** [(SRR6512652)](https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR6512652)

You can download and unzip input file with SRAtoolkit in `bio_tool/`:
```
cd examples
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR6512652/SRR6512652
../bio_tool/sratoolkit/bin/fastq-dump SRR6512652
cd ..
```

To execute the pipeline, please enter the following command:
```
sh run.sh examples/SRR6512652.fastq examples/piRNA_WS275.fa examples/mRNA_WS275.fa chira pirScan site
```
