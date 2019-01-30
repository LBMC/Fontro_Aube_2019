# Description

This program perform enrichment analyzes on frequencies of hexanucleotides, dinucleotides, nucleotides, nucleotides at a given codon position, codons, amino acids, and amino acid physicochemical features in exon sets.

The complete description on how the frequencies of the different features are computated is detailled in the sections entitled  "*Frequency of hexanucleotides, dinucleotides, nucleotides, codons, amino acids, and amino acid physicochemical features in exon sets*" and  "*Generation of sets of control exons and statistical analyses*"

# Prerequisites

## Installation of python and required modules

This program is written in **python 2.7** and was run in **Ubuntu 16.04**.

To install python 2.7 and python-pip you can run the following command (on ubuntu):

```sh
sudo apt install python2.7 python-pip
```

You can run the following command line to install the required modules :
```sh
pip install -r requirement.txt # installation of required python modules
```

## Building required files

### File `config.py`

Before doing anything, you must create a file named `config.py` with the following content

```py
# Connection to fasterDB
user="<user>"
password="<password>"
host="<host>"
database="<database>"
```
"\<user\>" and "\<password\>" corresponds to you username and password that allow you to connect the "\<host\>" with the FasterDB ("\<database\>") database on it.

### Frequencies files

To build the frequencies files, you must have created and filled the file `config.py` as described in the previous section.
Then just run the following command :

```sh
python2 frequences/frequencies_files_maker.py
mkdir nucleotides_riche_codon_ACE_CCE
python2 ACE_CCE_dic_maker.py --input frequences/input/input_ACE.txt --input2 frequences/input/input_CCE.txt --output nucleotides_riche_codon_ACE_CCE/
```

### Files for the summary figures.

If you want to use the option `--recap_figure yes` then you must run the following command :

```sh
mkdir frequence_recap_graphics
python2 getting_ace_cce_info.py frequence_recap_graphics
```

# Usage

To understand how the program works you can run the following command :

```sh
python main_project.py --help
```
