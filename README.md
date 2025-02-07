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
sh MutaCLASH.sh --input <input file> --regulator <regulator file> --transcript <transcript file> [--len <min hybrid length>] [--slen <max hybrid length>] [--link <adapter sequence>] [--trim <phred score>]
```
**required arguments:**
- **input file**: NGS data in FASTQ format.
- **regulator file**: regulator file in FASTA format.
- **transcript file**: transcript file in FASTA format.

**optional arguments(preprocessing):**
- **len**: Minimum hybrid length (default: 17).
- **slen**: Maximum hybrid length (default: 70).
- **link**: Adapter sequence (default: "None").
- **trim**: Phred score (default: 30).

MutaCLASH.sh will generate two .csv files: one for users to view the results of the MutaCLASH analysis (filename ending with _short), and the other for use with run_additional.sh to generate abundance results and figures.

_short.csv with columns `CLASH read sequence,read count,Deletion Sites on mRNA (Absolute Positions),Mismatch Sites on mRNA (Absolute Positions),Target RNA Name,Regulator RNA Name,Target RNA Region Found in CLASH Read,Region on CLASH Read identified as Target RNA,Region on CLASH Read identified as Regulator RNA,Regulator RNA Region Found in CLASH Read,pirScan score,Binding Energy Calculated by miRanda,Target Binding Sequence (miRanda),Regulator Binding Sequence (miRanda),RNAup Binding Energy,Target Binding Sequence (RNAup),Regulator Binding Sequence (RNAup)`

_short.csv example:
| CLASH read sequence | Read count | Deletion Sites on mRNA (Absolute Positions) | Mismatch Sites on mRNA (Absolute Positions) | Target RNA Name | Regulator RNA Name | Target RNA Region Found in CLASH Read | Region on CLASH Read identified as Target RNA | Region on CLASH Read identified as Regulator RNA | Regulator RNA Region Found in CLASH Read | pirScan score | Binding Energy Calculated by miRanda | Target Binding Sequence (miRanda) | Regulator Binding Sequence (miRanda) | RNAup Binding Energy | Target Binding Sequence (RNAup) | Regulator Binding Sequence (RNAup) |
|---------------------|------------|-----------------------------------------------|---------------------------------------------|----------------|----------------|-------------------------------------------|----------------------------------------------|----------------------------------------------|------------------------------------------|-------------|-------------------------------------|----------------------------------|----------------------------------|------------------|--------------------------------|--------------------------------|
| AAAAACACCGTCTTCCTCCAGTGGAGGCCTGGTTGTTTG | 6 | [] | [] | Y45F10D.12.1 | Y40H7A.12b | 427-446 | 2-21 | 1-18 | 22-39 | -27.5 | -13.83 | --gaAAAC-ACC-GTCTTCCt | cgtgTTTGTTGGTCCGGAGGt | -14.54 | --GAAAAC-ACCGTCTTCCTCCA | CGTGTTTGTTGGT--CCGGAGGT |
| AAAAACATCCATGCCCTCCAATCGTATTGGAGGCCTGGTTGTTTG | 3 | [] | [334] | C27A2.3.1 | Y40H7A.12b | 320-343 | 1-24 | 1-18 | 28-45 | -28.5 | -19.44 | --AAAAACATCCATGCTCTCCa | cgTGTTTGTTGGT-CCGGAGGt | -18.29 | --AAAAACATCCATGCTCTCCA | CGTGTTTGTTGGTCCG-GAGGT |


### run_additional.sh
```bash
sh run_additional.sh --input <PRG-1 or ALG-1>
```
**required arguments:**
- **input**: PRG-1 or ALG-1.


run_additional.sh is used to regenerate figures on paper; the user must first place ALG-1.csv and PRG-1.csv in the data/input folder before running this program.

[Reference metadata files and output files](http://nas.csblab.ee.ncku.edu.tw:32200/sharing/jSirL0jvo)

After executing the command, the pipeline will run and complete all the necessary steps. Please refer to the [examples](https://github.com/lu1215/MutaCLASH/tree/master/examples) we provided.


## Output
The output files are stored in the `data/output/` directory. The directory contains the following files:
- **CSV file**: Contains results with all information fields.
- **Figures:** The final generated figures are stored in the `figure/` subdirectory.
- **Logs**: Records commands, and summarizes the quantity, proportion, and distribution of various mutations, which are stored in the `log/` subdirectory.
- **Intermediate Files:** The intermediate files generated during the analysis are stored in various formats (.csv, etc.) and can be found in their respective tool directories.

### Figures
The output figures generated by the MutaCLASH pipeline include:

- **Score Distribution:** Presents the quantity and trend of different scores.
<img src="examples/fig/score.png" width=300 />

- **Mutation Distribution:** Provides information on the distribution of mutations, including deletions and substitutions.
<img src="examples/fig/distribution.png" width=300 />

- **Pairing Ratio:** Calculates and analyzes the pairing ratios at both global and individual coordinates.  
In statistical testing, ** and * indicate significant differences, with U-test P<0.05 and 0.10, respectively.
<img src="examples/fig/pairing_ratio_diff.png" width=600 />
<img src="examples/fig/pairing_ratio.png" width=600 />
<img src="examples/fig/pairing_ratio_at_position.png" width=300 />

- **Abundance Analysis:** Performs abundance analysis, comparing wild-type samples and fold-change measurements.
<!--img src="examples/fig/22G.png" width=300 /-->
<img src="examples/fig/fold_change.png" width=300 />
<img src="examples/fig/fold_change_per-score.png" width=500 />

- **Cumulative Distribution Function (CDF):** Calculates and visualizes the cumulative distribution function.
<!--img src="examples/fig/22G_CDF.png" width=300 /-->
<img src="examples/fig/fold_change_CDF.png" width=300 />

<!--Please refer to the corresponding tool documentation for more details on the specific output files and their interpretations.-->

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
