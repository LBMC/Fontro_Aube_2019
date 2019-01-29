# Description

The goal of this script is to make the wilcoxon-statistical test (against control exons) for the :
* Figure 1A : on SRSF1 regulated (repressed and activated) exons
* Figure 2B: on SRSF1 regulated exons (test for the enrichment of Guanine at codon position 1-2 and 3)
* Figure 3B: on SRSF1 regulated exons (test for the enrichment of Alanine, Glycine, Lysine and Proline encoded by SRSF1 regulated exons)
* Figure 5F: on hnRNPL and PTBP1 repressed exons (test for the enrichment in codon encoding hydroxyl-containing or negatively charged amino acids)

# Prerequisites

The script named `make_wilcoxon_text.py` uses **python 2.6** and the following modules:
* pandas
* rpy2
* subprocess
* numpy
* os

Make sure those are installed on your system.

Then you have to prepare the `data/` folder. You must create in this folder 3 sub-folders named `hnRNPL`, `PTBP1` and `SRSF1`.
For each of those folders, sub-folders (corresponding to the result of an enrichment analysis made by the tRNA program on a particular splicing factor in a particular cell line) with the following name must be created:
* N_n_SF_CELL_ID

Where:
* N and n are numbers (the value associated doesn't matter)
* SF is the name of the splicing factor studied : it must be the name of the parent folder
* CELL is a cell line
* ID is the database id from wich the data were downloaded before the enrichment analysis

Each subfolder much have the following content:
```sh
 N_n_SF_CELL_ID\
  |
  Exon_analysis_CCE/
    |
    down/
      |
      enrichment_report.xlsx
    up/
      |
      enrichment_report.xlsx
  ...
```
The folder "Exon_analysis_CCE/" is produced by the tRNA program with the following command line:

```sh
../tRNA_program/main_project.py --up input_repressed_exons_SF.txt --down input_activated_exon_SF.txt --exon_type CCE --enrichment True
```
`input_repressed_exons_SF.txt` and `input_activated_exon_SF.txt` are files containing repressed and activated exons by a splicing factor SF respectively.

# Launching the program

To perform the Mann-Withney Wilcoxon test, you only have to run:

```sh
python3 src/make_wilcoxon_text.py
```
