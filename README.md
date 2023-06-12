# Crosslink Induced Mutation Site (CIMS)

## Description
The **Crosslink Induced Mutation Site (CIMS)** project is designed to detect the coordinates of cims (crosslink induced mutation sites) in CLASH or iCLIP NGS data. It provides a comprehensive analysis pipeline for identifying mutation sites and binding sites in hybrid-reads derived from CLASH or iCLIP experiments.

## Features
- Utilizes CLASH Analyst to identify suitable NGS reads for analysis
- Uses ChiRA to identify suitable hybrid-reads
- Detects mutation information using Bowtie2
- Utilizes algorithms such as pirScan, miRanda, and RNAup to identify binding sites
- Generates visualizations of the distribution of mutations

## Usage
To run the CIMS pipeline, execute the following command:
```
sh run.sh <input_data.fa> <regulator.fa> <transcript.fa> <algorithm> <enrichment analysis type>
```
After executing the command, the pipeline will run and complete all the necessary steps.

## Output Files
The output files are stored in the `data/output/` directory. The directory contains the following files:
- Figures: The final generated figures are stored in the `figure/` subdirectory.
- Intermediate Files: The intermediate files generated during the analysis are stored in various formats (.csv, etc.) and can be found in their respective tool directories.

Please refer to the corresponding tool documentation for more details on the specific output files and their meanings.
