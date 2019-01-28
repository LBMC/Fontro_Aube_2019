#Prerequisites

Before doing anything you must create a file in your `scr/` a file named `config.py` with the following content :

```py
# Connection to fasterDB
user="<user>"
password="<password>"
host="<host>"
database="<database>"
```
"\<user\>" and "\<password\>" corresponds to you username and password that allow you to connect the "\<host\>" with the FasterDB ("\<database\>") database on it.

Then make sure that `python 2.7` is installed along with the following module:
* mysql.connector (version 2.4.1)

And run the following commands (this operation may take some time):

```sh
python2 src/make_control_dictionaries.py
```


# I - Program `fasta_generator.py`

## Description

Fasta generator program, (`fasta_generator.py`) allow to create a fasta file (containing DNA sequences) with your desired number of sequences.
 1. You can specify the proportion of nucleetides (alanine, guanine, thymine and cytosine) in your sequences.
 2. You can specify the proportion of dinucleotide (alanine, guanine, thymine and cytosine) in your sequences.
 3.  You can also create a sequences having codon proportion clause to the one of ACE/CCE/ALL exons of fasterDB and enriched it (**not impoverished**) with one particular di-nucleotide


## Prerequisites


This program works with `python 2.7` and the following modules :
* `random` : for random suffling of sequences
* `argparse` : for parsing arguments
* `os` : to check if a directory exists
* `math` : to truncate numbers
* `copy` : to make true copy of variables in python
* `sys` : to specify another directory of import


## Running the program

A description of all possible parameters can be displayed by running this command :

```sh
python2 src/fasta_generator.py --help
```

### Example of command

```sh
# Creation of a fasta file my_fasta.fasta in your current working directory
# It will contain 100 sequences having a size between 50 and 300 nt.
# Every sequences will have a proportion of alanine equals to 80%
# and a proportion to G, T, C equal to 6.67 %
python2 src/fasta_generator.py --prop_A 0.8 --filename "my_fasta" --nbr_seq
100 --size_inf 50 --size_max 300

# Creation of a fasta file my_fasta.fasta in your current working directory
# It will contain 100 sequences having a size between 50 and 300 nt.
# Every sequences will have a proportion of alanine equals to 50%
# and a proportion to G and  T equal to 20 %
# and a proportion of C of 10%
python2 src/fasta_generator.py --prop_A 0.5 --prop_C 0.1 --prop_T 0.2
--filename "my_fasta" --nbr_seq 100 --size_inf 50 --size_max 300

# Creation of a fasta file my_fasta.fasta in your current working directory
# It will contain 100 sequences having a size between 50 and 300 nt.
# Every sequences will have a proportion of AA equals to 50%
python2 src/fasta_generator.py --AA 0.5
--filename "my_fasta" --nbr_seq 100 --size_inf 50 --size_max 300

# Creation of a fasta file my_fasta.fasta in your current working directory
# It will contain 100 sequences having a size between 50 and 300 nt.
# Every sequences will have a proportion of AA equals to 50%.
# The sequence will be first generated with a codon prorotion equals to the one of
# CCE exons and then becom enriched in AA dinculeotide until it reach 50%

python2 src/fasta_generator.py --dnt AA 0.5 --ctrl CCE
--filename "my_fasta" --nbr_seq 100 --size_inf 50 --size_max 300

```

# II - Program `fasta_reverse_generator.py`

## Description

Fasta generator program, (`fasta_generator.py`) allow to create a fasta file (containing DNA sequences) with your desired number of sequences.
* You can specify the frequency of a feature and it will create a nucleotide sequence encoding for amino acid correponding to the feature choosen. Those amino acid will be present in the sequence according to the proportion choosen by the user.

## Prerequisites


This program works with `python 2.7` and the following modules :
* `random` : for random suffling of sequences
* `argparse` : for parsing arguments
* `os` : to check if a directory exists
* `sys` : to specify another directory of import
* `src/dictionary.py` : a file containing information about feature/amino acids.

You must launch the following comand before launching it to create the contol dictionaries
```sh
python src/make_control_dictionaries.py.py
```


## Running the program

A description of all possible parameters can be displayed by running this command :

```sh
python2 src/fasta_reverse_generator.py --help
```


### Example of command

```sh
# This command will cerate a fille in result/CCE_Hydroxilic_25.fasta. This file will contain 300 sequences having a length between 50 and 300 nt. Those sequence will encode 25% of hydroxilic amino acid with the same codon usage as the CCE exons of fasterDB
python src/fasta_reverse_generator.py --feature Hydroxylic --prob 25 --ctrl CCE --output result/

```

# III - Program make_comparison.py

## Description

The `make_comparison.py` programme allow to test if the sequences in a fasta file are enriched in some hexanucleotides, dinucleotides, codon, nucleotide at each codon position, amino acids and in different proteins features.
Hypergometric tests are made.

## Prerequisites


This program works with `python 2.7` and the following modules :
* `random` : for random suffling of sequences
* `Bio::SeqIO` : to parse fasta content
* `scipy.stats::hypergeom` : to make hypergeom tests
* `argparse` : for parsing arguments
* `os` : to check if a directory exists
* `sys` : to specify another directory of import
* `src/dictionary.py` : a file containing information about feature/amino acids.
* `re` : to allow the usage of regulard expression
* `rpy2` : to execute R functions
* `xlsxwriter` : to create an excel enrichment report
* `math` : to use the math pow function


You must launch the following comand before launching it to create the contol dictionaries
```sh
python src/make_control_dictionaries.py.py
```


## Running the program

A description of all possible parameters can be displayed by running this command :

```sh
python2 src/make_comparison.py --help
```

### Example of command

```sh

# This will compare a set of sequence in the fasta file with expeted value of a seque a CCE exons sequences and see if there is any enrichement in hexanucleotides/di-nucleotides/codon/nucleotides on each codon position/ in amino acids and in protein feature
python src/make_comparison.py --output result/ --ctrl CCE --fasta result/CCE_Hydroxylic_25.fasta --motif CT

```
