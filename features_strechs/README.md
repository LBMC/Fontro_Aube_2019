# Description

The goal of this script is to create the figures of the low complexity sequences in exons activated or represed by SRSF1, SRSF2, SRSF3 and TRA2.

# Prerequisites

The scripts `control_stretch_maker.py` and `stretch_evaluator.py` work with **python 3** and uses the following modules:
* pymysql
* os
* biopython
* pandas
* matplotlib
* argparse


## Control data files

>This part can be skipped if you have already a folder named control_dic into your `src` folder.

To produce the low complexity sequences charts you must first create a file name `config.py` having the following content:

```py
user = <username>
password = <password>
host = <host>
db = <fasterdb>
```
Those are informations that will allow you to connect the fasterDB database and thus retrieve the low complexity sequences for control exons.

Then you have to execute this command line :
```sh
python3 src/control_stretch_maker.py
```
This should create a folder name `src/control_dir` with control data files in it.

## preparing your `data` folder

Your data folder must contained a folder named `SF_regulated_exons` having the following structure :
```sh
SF_regulated_exon/
  |
  folder1
    query_resuls_up_exons.xlsx
    query_resuls_down_exons.xlsx
```
The `query_resuls_up_exons.xlsx` and `query_resuls_down_exons.xlsx` must have been created with the tRNA program.
The `query_resuls_down_exons.xlsx` should contains only exon that are activated by a splicing factor  (i.e the tRNA program that produced this file was launched with an input containing exon activated by a splicing factor)
The `query_resuls_up_exons.xlsx` should contains only exon that are repressed by a splicing factor  (i.e the tRNA program that produced this file was launched with an input containing exon repressed by a splicing factor)

To produce those files with the tRNA program you can launch :

```sh
python2 main_project.py --up input_up.txt --down input_down.txt --enrichment False
```

# Executing the program

you can launch the following command:

```sh
bash src/command.sh
```

It will produce some figures that have been used in the article.
