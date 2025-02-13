# MutaCLASH Data Processing Guide

## Before Analysis
To get started, clone the repository and install the required dependencies:

```bash
git clone https://github.com/lu1215/MutaCLASH.git
cd MutaCLASH/
apt-get install -y samtools bowtie2
pip install -r requirements.txt
```

## PRG-1 CLASH Data Processing
1. **[Download the PRG-1.csv file](http://nas.csblab.ee.ncku.edu.tw:32200/fsdownload/jSirL0jvo/example_data_for_github)** and place it in the following directory:
   ```
   MutaCLASH/data/input
   ```
2. Run the processing script:
   ```bash
   sh run_additional.sh --input PRG-1
   ```
3. The output will be saved in:
   ```
   MutaCLASH/data/output/PRG-1_<date_time>/
   ```
   This folder contains:
   - Two CSV files:
     - A detailed dataset with MutaCLASH and abundance information.
     - A `_short.csv` file with simplified, more interpretable columns.
   - A `figure` folder containing:
     - Four subfolders: `abu_plot`, `distribution_plot`, `G22_plot`, and `pairing_ratio_plot`.
     - A score plot.
   
   For a detailed explanation of each output file and figure, refer to the [MutaCLASH GitHub documentation](https://github.com/lu1215/MutaCLASH/tree/master?tab=readme-ov-file#output).

---

## ALG-1 CLASH Data Processing
1. **[Download the ALG-1.csv file](http://nas.csblab.ee.ncku.edu.tw:32200/fsdownload/jSirL0jvo/example_data_for_github)** and place it in:
   ```
   MutaCLASH/data/input
   ```
2. Execute the following command:
   ```bash
   sh run_additional.sh --input ALG-1
   ```
3. The results will be saved in:
   ```
   MutaCLASH/data/output/ALG-1_<date_time>/
   ```
   This folder contains:
   - Two CSV files:
     - A comprehensive dataset with MutaCLASH and abundance data.
     - A `_short.csv` version with simplified column names for easier interpretation.
   - A `figure` directory, which includes:
     - Subdirectories: `abu_plot`, `distribution_plot`, `G22_plot`, and `pairing_ratio_plot`.
     - A score plot.

For further details on the figures and output structure, visit the [MutaCLASH GitHub documentation](https://github.com/lu1215/MutaCLASH/tree/master?tab=readme-ov-file#output).

