# MutaCLASH

## Description
The **MutaCLASH** project is designed to detect the coordinates of Crosslink Induced Mutation Sites (CIMS) in NGS data. It provides a comprehensive analysis pipeline for identifying mutation sites and binding sites in hybrid-reads derived from CLASH or iCLIP experiments.

## Features
- Uses [ChiRA](https://github.com/pavanvidem/chira) to identify suitable hybrid-reads.
- Detects mutation information using [Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml) and [BWA](https://bio-bwa.sourceforge.net).
- Utilizes algorithms such as [pirScan](http://cosbi4.ee.ncku.edu.tw/pirScan/), [miRanda](https://bioweb.pasteur.fr/packages/pack@miRanda@3.3a), and [RNAup](https://github.com/ViennaRNA/ViennaRNA) to identify binding sites.
- Generates visualizations of the distribution of mutations.

## Usage
### MutaCLASH.sh
To run only MutaCLASH pipeline, execute the following command:
```bash
sh MutaCLASH.sh --input <input file> --regulator <regulator file> --transcript <transcript file> [--len <min hybrid length>] [--slen <max hybrid length>] [--link <adapter sequence>] [--trim <phred score>] [--transposon]
```
**required arguments:**
- **input file**: NGS data in FASTQ format. (relative address of MutaCLASH folder)
- **regulator file**: regulator file in FASTA format. (relative address of MutaCLASH folder)
- **transcript file**: transcript file in FASTA format. (relative address of MutaCLASH folder)

**optional arguments(preprocessing):**
- **len**: Minimum hybrid length (default: 17).
- **slen**: Maximum hybrid length (default: 70).
- **link**: Adapter sequence (default: "None").
- **trim**: Phred score (default: 30).
- **transposon**: Setting BWA-MEM in ChiRA with `-a` parameter, to retain equally scoring alignments (-a in BWA-MEM), this parameter can be used to analyse transposon RNA data.

MutaCLASH.sh will generate two .csv files: one for users to view the results of the MutaCLASH analysis (filename ending with _short), and the other for use with run_additional.sh to generate abundance results and figures.

_short.csv with columns `CLASH read sequence,read count,Target RNA Name,Regulator RNA Name,Target RNA Region Found in CLASH Read,Region on CLASH Read identified as Target RNA,Region on CLASH Read identified as Regulator RNA,Regulator RNA Region Found in CLASH Read,Deletion Sites on mRNA (Absolute Positions),Mismatch Sites on mRNA (Absolute Positions),pirScan score,pirScan binding site,pirScan Target RNA sequence,pirScan Regulator RNA sequence,Binding Energy Calculated by miRanda,miRanda binding site,Target Binding Sequence (miRanda),Regulator Binding Sequence (miRanda),RNAup Binding Energy,RNAup binding site,Target Binding Sequence (RNAup),Regulator Binding Sequence (RNAup)`

_short.csv example:
<!-- | CLASH read sequence                                   | read count | Target RNA Name | Regulator RNA Name | Target RNA Region Found in CLASH Read | Region on CLASH Read identified as Target RNA | Region on CLASH Read identified as Regulator RNA | Regulator RNA Region Found in CLASH Read | Deletion Sites on mRNA (Absolute Positions) | Mismatch Sites on mRNA (Absolute Positions) | pirScan score | pirScan binding site | pirScan Target RNA sequence | pirScan Regulator RNA sequence | Binding Energy Calculated by miRanda | miRanda binding site | Target Binding Sequence (miRanda) | Regulator Binding Sequence (miRanda) | RNAup Binding Energy | RNAup binding site | Target Binding Sequence (RNAup) | Regulator Binding Sequence (RNAup) |
|-------------------------------------------------------|------------|-----------------|---------------------|-----------------------------------------|-----------------------------------------------|-----------------------------------------------|---------------------------------------------|-----------------------------------------------|-----------------------------------------------|--------------|----------------------|-----------------------------|-----------------------------|---------------------------------|--------------------|----------------------------------|----------------------------------|------------------|----------------|--------------------------------|--------------------------------|
| AAAAACACCGTCTTCCTCCAGTGGAGGCCTGGTTGTTTG              | 6          | Y45F10D.12.1    | Y40H7A.12b          | 427-446                                 | 2-21                                         | 1-18                                         | 22-39                                       | []                                            | []                                            | -27.5        | 425-445              | GACCTCCTTCTGCCACAAAAG       | TGGAGGCCTGGTTGTTTGTGC       | -13.83                          | 425-446            | &nbsp;--gaAAAC-ACC-GTCTTCCt         | cgtgTTTGTTGGTCCGGAGGt         | -14.54            | 425-446        | &nbsp;--GAAAAC-ACCGTCTTCCTCCA     | CGTGTTTGTTGGT--CCGGAGGT     |
| AAAAACATCCATGCCCTCCAATCGTATTGGAGGCCTGGTTGTTTG        | 3          | C27A2.3.1       | Y40H7A.12b          | 320-343                                 | 1-24                                         | 1-18                                         | 28-45                                       | []                                            | [334]                                         | -28.5        | 321-341              | CTAACCTCTCGTACCTACAAA       | TGGAGGCCTGGTTGTTTGTGC       | -19.44                          | 319-342            | &nbsp;--<br>AAAAACATCCATGCTCTCCa       | cgTGTTTGTTGGT-CCGGAGGt       | -18.29            | 319-342        | &nbsp;--<br>AAAAACATCCATGCTCTCCA     | CGTGTTTGTTGGTCCG-GAGGT     |
| AAAAACATCCATGCTCTCCAATCGACACTGCAAACTATTGAGGCCTGGTTGTTTG | 2          | C27A2.3.1       | Y40H7A.12b          | 320-354                                 | 1-35                                         | 3-18                                         | 40-55                                       | []                                            | []                                            | -28.5        | 331-351              | AACGTCACAGCTAACCTCTCG       | TGGAGGCCTGGTTGTTTGTGC       | -19.44                          | 319-353            | &nbsp;--<br>AAAAACATCCATGCTCTCCa       | cgTGTTTGTTGGT-CCGGAGGt       | -16.89            | 319-353        | &nbsp;--<br>AAAAACATCCATGCTCTCCA     | CGTGTTTGTTGGTCCG-GAGGT     | -->


<table border="1">
<tr>
<th>CLASH read sequence</th>
<th>read count</th>
<th>Target RNA Name</th>
<th>Regulator RNA Name</th>
<th>Target RNA Region Found in CLASH Read</th>
<th>Region on CLASH Read identified as Target RNA</th>
<th>Region on CLASH Read identified as Regulator RNA</th>
<th>Regulator RNA Region Found in CLASH Read</th>
<th>Deletion Sites on mRNA (Absolute Positions)</th>
<th>Mismatch Sites on mRNA (Absolute Positions)</th>
<th>pirScan score</th>
<th>pirScan binding site</th>
<th>pirScan Target RNA sequence</th>
<th>pirScan Regulator RNA sequence</th>
<th>Binding Energy Calculated by miRanda</th>
<th>miRanda binding site</th>
<th>Target Binding Sequence (miRanda)</th>
<th>Regulator Binding Sequence (miRanda)</th>
<th>RNAup Binding Energy</th>
<th>RNAup binding site</th>
<th>Target Binding Sequence (RNAup)</th>
<th>Regulator Binding Sequence (RNAup)</th>
</tr>
<tr>
<td nowrap="nowrap" >AAAAACACCGTCTTCCTCCAGTGGAGGCCTGGTTGTTTG</td>
<td nowrap="nowrap" >6</td>
<td nowrap="nowrap" >Y45F10D.12.1</td>
<td nowrap="nowrap" >Y40H7A.12b</td>
<td nowrap="nowrap" >427-446</td>
<td nowrap="nowrap" >2-21</td>
<td nowrap="nowrap" >1-18</td>
<td nowrap="nowrap" >22-39</td>
<td nowrap="nowrap" >[]</td>
<td nowrap="nowrap" >[]</td>
<td nowrap="nowrap" >-27.5</td>
<td nowrap="nowrap" >425-445</td>
<td nowrap="nowrap" >GACCTCCTTCTGCCACAAAAG</td>
<td nowrap="nowrap" >TGGAGGCCTGGTTGTTTGTGC</td>
<td nowrap="nowrap" >-13.83</td>
<td nowrap="nowrap" >425-446</td>
<td nowrap="nowrap" >--gaAAAC-ACC-GTCTTCCt</td>
<td nowrap="nowrap" >cgtgTTTGTTGGTCCGGAGGt</td>
<td nowrap="nowrap" >-14.54</td>
<td nowrap="nowrap" >425-446</td>
<td nowrap="nowrap" >--GAAAAC-ACCGTCTTCCTCCA</td>
<td nowrap="nowrap"  >CGTGTTTGTTGGT--CCGGAGGT</td>
</tr>
<tr>
<td nowrap="nowrap" >AAAAACATCCATGCCCTCCAATCGTATTGGAGGCCTGGTTGTTTG</td>
<td nowrap="nowrap" >3</td>
<td nowrap="nowrap" >C27A2.3.1</td>
<td nowrap="nowrap" >Y40H7A.12b</td>
<td nowrap="nowrap" >320-343</td>
<td nowrap="nowrap" >1-24</td>
<td nowrap="nowrap" >1-18</td>
<td nowrap="nowrap" >28-45</td>
<td nowrap="nowrap" >[]</td>
<td nowrap="nowrap" >[334]</td>
<td nowrap="nowrap" >-28.5</td>
<td nowrap="nowrap" >321-341</td>
<td nowrap="nowrap" >CTAACCTCTCGTACCTACAAA</td>
<td nowrap="nowrap" >TGGAGGCCTGGTTGTTTGTGC</td>
<td nowrap="nowrap" >-19.44</td>
<td nowrap="nowrap" >319-342</td>
<td nowrap="nowrap" >--AAAAACATCCATGCTCTCCa</td>
<td nowrap="nowrap" >cgTGTTTGTTGGT-CCGGAGGt</td>
<td nowrap="nowrap" >-18.29</td>
<td nowrap="nowrap" >319-342</td>
<td nowrap="nowrap" >--AAAAACATCCATGCTCTCCA</td>
<td nowrap="nowrap" >CGTGTTTGTTGGTCCG-GAGGT</td>
</tr>
<tr>
<td nowrap="nowrap" >AAAAACATCCATGCTCTCCAATCGACACTGCAAACTATTGAGGCCTGGTTGTTTG</td>
<td nowrap="nowrap" >2</td>
<td nowrap="nowrap" >C27A2.3.1</td>
<td nowrap="nowrap" >Y40H7A.12b</td>
<td nowrap="nowrap" >320-354</td>
<td nowrap="nowrap" >1-35</td>
<td nowrap="nowrap" >3-18</td>
<td nowrap="nowrap" >40-55</td>
<td nowrap="nowrap" >[]</td>
<td nowrap="nowrap" >[]</td>
<td nowrap="nowrap" >-28.5</td>
<td nowrap="nowrap" >331-351</td>
<td nowrap="nowrap" >AACGTCACAGCTAACCTCTCG</td>
<td nowrap="nowrap" >TGGAGGCCTGGTTGTTTGTGC</td>
<td nowrap="nowrap" >-19.44</td>
<td nowrap="nowrap" >319-353</td>
<td nowrap="nowrap" >--AAAAACATCCATGCTCTCCa</td>
<td nowrap="nowrap" >cgTGTTTGTTGGT-CCGGAGGt</td>
<td nowrap="nowrap" >-16.89</td>
<td nowrap="nowrap" >319-353</td>
<td nowrap="nowrap" >--AAAAACATCCATGCTCTCCA</td>
<td nowrap="nowrap" >CGTGTTTGTTGGTCCG-GAGGT</td>
</tr>
</table>

## Requirements
Running MutaCLASH require Linux or MacOS. Other Unix environments will probably work but have not been tested. Windows users can use [Windows Subsystem for Linux](https://docs.microsoft.com/en-us/windows/wsl/install-win10).

To install some necessary tools and packages, execute the following command:
```bash
$ apt-get install -y samtools bowtie2
$ pip install -r requirements.txt
```

- SAMtools >= 0.1.19
- Bowtie2 >= 2.4.0
- Python >= 3.5
- bcbio-gff >= 0.6.9
- biopython >= 1.76
- cutadapt >= 2.10
- matplotlib >= 2.2.2
- numpy >= 1.12.1
- pandas >= 0.23.0
- pysam >= 0.20.0
- scipy >= 1.1.0
- seaborn >= 0.9
- statannot = 0.2.3
- tqdm >= 4.64.0
- xlrd >= 1.2.0

## Docker
If you have any concerns about environment setup, feel free to use the Docker version directly.
```bash
$ docker pull ryanccj/mutaclash
$ docker run -it ryanccj/mutaclash
```

Or you can choose to built from this project.
```bash
$ docker build -t <image_name> .
$ docker run -it <image_name>
```

## LICENSE
Please refer to our [MIT license](https://github.com/RyanCCJ/MutaCLASH/blob/master/LICENSE).
