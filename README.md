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
sh MutaCLASH.sh --input <input file> --regulator <regulator file> --transcript <transcript file> [--prep <preprocessing tool>] [--len <min hybrid length>] [--slen <max hybrid length>] [--link <adapter sequence>] [--trim <phred score>] [--transposon]
```
**required arguments:**
- **input file**: NGS data in FASTQ format. (relative address of MutaCLASH folder)
- **regulator file**: regulator file in FASTA format. (relative address of MutaCLASH folder)
- **transcript file**: transcript file in FASTA format. (relative address of MutaCLASH folder)

**optional arguments:**
- **about preprocessing**
    - **prep**: Prepprocessing tool, User can choose "cutadapt" or "trimgalore" (default: trimgalore), if choose cutadapt, `--link` will become required argument.
    - **len**: Minimum hybrid length (default: 17).
    - **slen**: Maximum hybrid length (default: 70).
    - **link**: Adapter sequence (default: "None").
    - **trim**: Phred score (default: 30).
- **transposon**: Setting BWA-MEM in ChiRA with `-a` parameter, to retain equally scoring alignments (-a in BWA-MEM), this parameter can be used to analyse transposon RNA data.

example command for analysing Worm PRG-1 data:
before executing command below, please put `PRG-1_rep1.fastq` into `MutaCLASH/data/input`.

```bash
sh MutaCLASH.sh --input PRG-1_rep1.fastq --regulator data/reference/piRNA_WS275_ori_order.fa --transcript data/reference/mRNA_WS275.fa --prep cutadapt --link AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
```

MutaCLASH.sh will generate two .csv files: one for users to view the results of the MutaCLASH analysis (filename ending with _short), and the other for use with run_additional.sh to generate abundance results and figures.

_short.csv with columns `CLASH read sequence,read count,normalized read count (evenly distributed),Target RNA Name,Regulator RNA Name,Target RNA Region Found in CLASH Read,Region on CLASH Read identified as Target RNA,Regulator RNA Region Found in CLASH Read,Region on CLASH Read identified as Regulator RNAd,Deletion Sites on mRNA (Absolute Positions),Mismatch Sites on mRNA (Absolute Positions),pirScan score,pirScan binding site,pirscan Target RNA sequence,pirscan Regulator RNA sequence,Binding Energy Calculated by miRanda,miRanda binding site,Target Binding Sequence (miRanda),Regulator Binding Sequence (miRanda),RNAup Binding Energy,RNAup binding site,Target Binding Sequence (RNAup),Regulator Binding Sequence (RNAup),single_read_mut_stat_significance_D,single_read_mut_stat_significance_M`

* `single_read_mut_stat_significance_D`, `single_read_mut_stat_significance_M`:
At the single-read level, these columns indicate whether the number of CIMS within the identified target region of this read is statistically significantly different from what would be expected based on the entire target RNA.
…_D evaluates deletions, and …_M evaluates substitutions (mismatches). Values are Boolean (True/False) for significance.

_short.csv example:

<table border="1">
<tr>
<th>CLASH read sequence</th>
<th>read count</th>
<th>normalized read count (evenly distributed)</th>
<th>Target RNA Name</th>
<th>Regulator RNA Name</th>
<th>Target RNA Region Found in CLASH Read</th>
<th>Region on CLASH Read identified as Target RNA</th>
<th>Regulator RNA Region Found in CLASH Read</th>
<th>Region on CLASH Read identified as Regulator RNAd</th>
<th>Deletion Sites on mRNA (Absolute Positions)</th>
<th>Mismatch Sites on mRNA (Absolute Positions)</th>
<th>pirScan score</th>
<th>pirScan binding site</th>
<th>pirscan Target RNA sequence</th>
<th>pirscan Regulator RNA sequence</th>
<th>Binding Energy Calculated by miRanda</th>
<th>miRanda binding site</th>
<th>Target Binding Sequence (miRanda)</th>
<th>Regulator Binding Sequence (miRanda)</th>
<th>RNAup Binding Energy</th>
<th>RNAup binding site</th>
<th>Target Binding Sequence (RNAup)</th>
<th>Regulator Binding Sequence (RNAup)</th>
<th>single_read_mut_stat_significance_D</th>
<th>single_read_mut_stat_significance_M</th>
</tr>
<tr>
<td nowrap="nowrap">AAAAAAAAAAAGAAAGATTTGTTGAAAGTTTCAACAATCTAATCATTTTA</td>
<td nowrap="nowrap">1</td>
<td nowrap="nowrap">1.0</td>
<td nowrap="nowrap">T10G3.5a.1</td>
<td nowrap="nowrap">21ur-9428</td>
<td nowrap="nowrap">1766-1782</td>
<td nowrap="nowrap">11-27</td>
<td nowrap="nowrap">1-21</td>
<td nowrap="nowrap">30-50</td>
<td nowrap="nowrap">[]</td>
<td nowrap="nowrap">[]</td>
<td nowrap="nowrap">-29.5</td>
<td nowrap="nowrap">1765-1785</td>
<td nowrap="nowrap">AGGTAAAGTTGTTTAGAAAGA</td>
<td nowrap="nowrap">TTCAACAATCTAATCATTTTA</td>
<td nowrap="nowrap">False</td>
<td nowrap="nowrap">False-False</td>
<td nowrap="nowrap">False</td>
<td nowrap="nowrap">False</td>
<td nowrap="nowrap">-10.71</td>
<td nowrap="nowrap">1761-1785</td>
<td nowrap="nowrap">---AAGGAGAAAGATTTGTTGAA</td>
<td nowrap="nowrap">ATTTTACTAATCTAACAACTT</td>
<td nowrap="nowrap">False</td>
<td nowrap="nowrap">False</td>
</tr>
<tr>
<td nowrap="nowrap">AAATGGAAGAGGAACGACAAAAGTCTTTTCGTTCCTCTATCTAAA</td>
<td nowrap="nowrap">2</td>
<td nowrap="nowrap">2.0</td>
<td nowrap="nowrap">T10G3.5a.1</td>
<td nowrap="nowrap">21ur-638</td>
<td nowrap="nowrap">1792-1814</td>
<td nowrap="nowrap">1-23</td>
<td nowrap="nowrap">1-21</td>
<td nowrap="nowrap">24-44</td>
<td nowrap="nowrap">[]</td>
<td nowrap="nowrap">[]</td>
<td nowrap="nowrap">-13.5</td>
<td nowrap="nowrap">1792-1812</td>
<td nowrap="nowrap">AAAACAGCAAGGAGAAGGTAA</td>
<td nowrap="nowrap">TCTTTTCGTTCCTCTATCTAA</td>
<td nowrap="nowrap">False</td>
<td nowrap="nowrap">False-False</td>
<td nowrap="nowrap">False</td>
<td nowrap="nowrap">False</td>
<td nowrap="nowrap">-19.46</td>
<td nowrap="nowrap">1791-1813</td>
<td nowrap="nowrap">ATGGAAGAGGAACGACAAAAG-</td>
<td nowrap="nowrap">AATCTATCTCCTTGCTTTTCT</td>
<td nowrap="nowrap">False</td>
<td nowrap="nowrap">False</td>
</tr>
<tr>
<td nowrap="nowrap">AAATTGAGTCTCTGAAAACTACGTTGTTTTCAGAAGCACAATTTA</td>
<td nowrap="nowrap">2</td>
<td nowrap="nowrap">2.0</td>
<td nowrap="nowrap">T10G3.5a.1</td>
<td nowrap="nowrap">21ur-12811</td>
<td nowrap="nowrap">2689-2711</td>
<td nowrap="nowrap">1-23</td>
<td nowrap="nowrap">1-21</td>
<td nowrap="nowrap">24-44</td>
<td nowrap="nowrap">[]</td>
<td nowrap="nowrap">[]</td>
<td nowrap="nowrap">-2.5</td>
<td nowrap="nowrap">2688-2708</td>
<td nowrap="nowrap">ATCAAAAGTCTCTGAGTTAAA</td>
<td nowrap="nowrap">TTGTTTTCAGAAGCACAATTT</td>
<td nowrap="nowrap">False</td>
<td nowrap="nowrap">False-False</td>
<td nowrap="nowrap">False</td>
<td nowrap="nowrap">False</td>
<td nowrap="nowrap">-14.26</td>
<td nowrap="nowrap">2688-2710</td>
<td nowrap="nowrap">AAATTGAGTCTCTGAAAACTA</td>
<td nowrap="nowrap">TTTAACACGAAGACTTTTGTT</td>
<td nowrap="nowrap">False</td>
<td nowrap="nowrap">False</td>
</tr>
</table>


## Analysis Workflow

The **MutaCLASH** pipeline is composed of multiple well-structured steps, integrating several bioinformatics tools to process raw NGS reads and identify mutation sites and binding interactions between RNAs. Below is an overview of the full analysis workflow:

1. **Preprocessing (Step 1)**  
   Raw reads in FASTQ format are first processed to remove low-quality bases and adapter sequences.  
   - By default, the pipeline uses **Trim Galore** for trimming.  
   - Alternatively, **Cutadapt** can be selected via the `--prep cutadapt` option (in this case, an adapter sequence must be provided via `--link`).  
   - The preprocessing step also includes **deduplication** to merge identical reads and to calculate readcount.

2. **Hybrid Read Detection using ChiRA (Step 2)**  
   Preprocessed reads are analyzed with [**ChiRA**](https://github.com/pavanvidem/chira), which identifies **chimeric reads** (hybrid reads) representing interactions between regulator and target RNAs.  
   - This step aligns reads to both regulator and transcript reference sequences.  
   - ChiRA uses **BWA-MEM** internally; enabling `--transposon` adds the `-a` flag to retain multiple equally scoring alignments, useful for repetitive/transposon data.

3. **Mutation (Deletion and Substitution) Identification (Step 3)**  
   Using the **find_deletion** module, the pipeline identifies **Crosslink Induced Mutation Sites (CIMS)** such as deletions and substitutions in target RNA regions from aligned hybrid reads.  
   - Detected mutation sites are recorded with their **absolute positions** in transcripts.

4. **Binding Site Prediction (Step 4)**  
   The detected hybrid reads and their mapped coordinates are further analyzed using three independent algorithms to predict potential **regulator–target binding sites**:  
   - [**pirScan**](http://cosbi4.ee.ncku.edu.tw/pirScan/): Predicts piRNA target sites based on conserved seed rules.  
   - [**miRanda**](https://bioweb.pasteur.fr/packages/pack@miRanda@3.3a): Calculates complementarity and binding energy between regulator and target.  
   - [**RNAup**](https://github.com/ViennaRNA/ViennaRNA): Estimates RNA–RNA binding energy considering secondary structure accessibility.  

5. **Single-Read Statistics (Step 5 — non-chimeric reads)**  
   In this stage, single-read refers to any read that ChiRA does not classify as chimeric (e.g., only one arm aligns or no fusion junction is detected).
   This step computes and exports CIMS (crosslink-induced mutation sites) from single-read data in the ChiRA output. These results serve as inputs for downstream comparisons and visualization.

6. **Data Processing and Integration (Step 6)**  
   - The results from mutation detection, binding predictions, and statistics are merged into a **comprehensive summary table**, ensuring consistent formatting and column naming.
   - If the same hybrid read appears in the dataset both with CIMS (deletions or substitutions) and without CIMS, we remove the row with CIMS and keep the CIMS-free row as the single canonical record for that read. This ensures each hybrid read is represented by its best (error-free) instance.

7. **Result Compilation (Step 7)**  
   All output files, logs, and configuration details are collected into a timestamped directory under `data/output/`.  
   - The main result is a `.csv` file with full annotations.  
   - A corresponding `_short.csv` version is also generated, containing key columns for quick inspection or downstream abundance visualization.  
   - Column names are reformatted for clarity (e.g., `read_count` → `read count`, `mir_score` → `miRanda score`).

At the end of the pipeline, you will obtain two main CSV files:
- `<DATA>.csv`: Full detailed output for downstream analysis and visualization.
- `<DATA>_short.csv`: Simplified summary.

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
