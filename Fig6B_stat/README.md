# Description

The goal of this script is perform the statistical analysis of the figure 6B of the article.

# Prerequisites.

The list of all activated and repressed exons by SRSF1, SRSF2, SRSF3, and TRA2 must have been processed by the tRNA program with this command line :
```sh
python3 main_program.py --up all_repressed_exon_SFname.txt --down all_activated_exon_SFname.txt --enrichment False --exon_type CCE
```

This command line will produce some query result_files that must be placed into the ``data`` folder of this project with the following architecture:
```sh
data/
  |
  N_ALL_SFname/
        |
        query_result_activated.xlsx
        query_result_repressed.xlsx
  ...
```

With ``N_ALL_SFname`` correponds to a folder in ``data``, where:
* N is a number
* SFname is the name of the splicing factor that regulate the exons contained in the files  ``query_result_activated.xlsx`` and ``query_result_repressed`` produced by the tRNA program.
* ALL mean that the list of exon that have been processed by the tRNA program corresponds to each exons activated or repressend by the splicing factor ``SFname`` in at leat one cell line. If an exon is regulated by ``SFname`` in multiple cell lines, it is keeped in the list if it is always regulated the same way.

The script works with ``python3`` and uses the following modules:
* pandas
* subprocess
* rpy2


# Execution

```sh
python3 src/Fig6B_stats_calculator.py
```

This will produce a file containing the p-values of the figure 6B in the `result` folder.
