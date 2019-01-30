####################################################
#                   Description
####################################################
# Programme that performs the creation of a file named result_ACE_CCE_boxplot_recap.py
# It will contains foreach codon and amino_acid the list of their proportion in ACE and CCE sets.
# Those list will contain as much values as the number of ACE and CC exons presents in fasterDB

####################################################
#                    Imports
####################################################

from dictionnary import dic_sequence
import mysql.connector
from exon_set_class import *
from list_exon import *
from dictionnary import amino_acid2nature
from dictionnary import codon2aminoAcid
import sys
import config

#####################################################
#                  Functions
#####################################################


def connection():
    """
    :return: an object that contains all the information you need to connect to fasterDB
    """
    cnx = mysql.connector.connect(user=config.user, password=config.password, host=config.host, database=config.database,
                                  buffered=True)
    return cnx


def print_progress(iteration, total, prefix='', suffix='', decimals=1, bar_length=100):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        barLength   - Optional  : character length of bar (Int)
    """
    format_str = "{0:." + str(decimals) + "f}"
    percent = format_str.format(100 * (iteration / float(total)))
    filled_length = int(round(bar_length * iteration / float(total)))
    bar = "#" * filled_length + '-' * (bar_length - filled_length)
    sys.stdout.write('\r%s |%s| %s%s %s' % (prefix, bar, percent, '%', suffix)),
    if iteration == total:
        sys.stdout.write('\n')
    sys.stdout.flush()


def creating_list_exon(exon_type):
    """
    :param exon_type: (string) the type of the exon we want to analyse
    :return: (list of exon instance)
    """
    cnx = connection()
    dic = dic_sequence(cnx, exon_type)
    set_result = ListExon()
    counter = 1
    total = len(dic.keys())
    pre_info = "Creating the " + str(exon_type) + " list..."
    for key in dic.keys():
        print_progress(counter, total, prefix=pre_info, suffix='', decimals=1, bar_length=100)
        exon_info = key.split("_")
        set_result.exon_list.append(ExonSetClass(exon_info))
        set_result.exon_list[-1].cds_sequence = dic[key][0]
        set_result.exon_list[-1].offset = dic[key][1]
        set_result.exon_list[-1].found_codon_amino_acid_and_aa_nature()
        counter += 1
    return set_result


def writing_results(output_folder):
    """
    :param output_folder: (string) path were the results will be created.

    1 - Write a file named result_ACE_CCE_boxplot_recap.py
    It will contains foreach codon and amino_acid the list of their proportion in ACE and CCE sets.
    Those list will contain as much values as the number of ACE and CC exons presents in fasterDB
    2 - Write a file named result_ACE_CCE_barplot_recap.py. It will contains foreach codon and amino_acid the proportion
    of the ACE and CCE exon that have 0, 1, 2, 3, 4, 5 or more of those codon, amino acid in their sequence
    """
    list_exon = creating_list_exon("ACE")
    my_file = open(output_folder + "/result_ACE_CCE_boxplot_recap.py", "w")
    barplot_file = open(output_folder + "/result_ACE_CCE_barplot_recap.py", "w")
    counter = 1
    total = len(codon2aminoAcid.keys())
    pre_info = "Creating the ACE codon proportion list..."
    for codon in codon2aminoAcid.keys():
        print_progress(counter, total, prefix=pre_info, suffix='', decimals=1, bar_length=100)
        res = list_exon.get_exon_proportion_of_codon(codon)
        res_count = list_exon.get_exon_proportion_count_of_codon(codon)
        my_file.write("ACE_" + str(codon) + " = " + str(res) + "\n")
        barplot_file.write("ACE_" + str(codon) + " = " + str(res_count) + "\n")
        counter += 1
    counter = 1
    total = len(amino_acid2nature.keys())
    pre_info = "Creating the ACE aa proportion list..."
    for aa in amino_acid2nature.keys():
        print_progress(counter, total, prefix=pre_info, suffix='', decimals=1, bar_length=100)
        res = list_exon.get_exon_proportion_of_amino_acid(aa)
        res_count = list_exon.get_exon_proportion_count_of_amino_acid(aa)
        if aa != "*":
            my_file.write("ACE_" + str(aa) + " = " + str(res) + "\n")
            barplot_file.write("ACE_" + str(aa) + " = " + str(res_count) + "\n")
        else:
            my_file.write("ACE_stop  = " + str(res) + "\n")
            barplot_file.write("ACE_stop  = " + str(res_count) + "\n")
        counter += 1
    list_exon = creating_list_exon("CCE")
    counter = 1
    total = len(codon2aminoAcid.keys())
    pre_info = "Creating the CCE codon proportion list..."
    for codon in codon2aminoAcid.keys():
        print_progress(counter, total, prefix=pre_info, suffix='', decimals=1, bar_length=100)
        res = list_exon.get_exon_proportion_of_codon(codon)
        res_count = list_exon.get_exon_proportion_count_of_codon(codon)
        my_file.write("CCE_" + str(codon) + " = " + str(res) + "\n")
        barplot_file.write("CCE_" + str(codon) + " = " + str(res_count) + "\n")
        counter += 1
    counter = 1
    total = len(amino_acid2nature.keys())
    pre_info = "Creating the CCE aa proportion list..."
    for aa in amino_acid2nature.keys():
        print_progress(counter, total, prefix=pre_info, suffix='', decimals=1, bar_length=100)
        res = list_exon.get_exon_proportion_of_amino_acid(aa)
        res_count = list_exon.get_exon_proportion_count_of_amino_acid(aa)
        if aa != "*":
            my_file.write("CCE_" + str(aa) + " = " + str(res) + "\n")
            barplot_file.write("CCE_" + str(aa) + " = " + str(res_count) + "\n")
        else:
            my_file.write("CCE_stop" + " = " + str(res) + "\n")
            barplot_file.write("CCE_stop  = " + str(res_count) + "\n")
        counter += 1
    my_file.close()

if __name__ == "__main__":
    if len(sys.argv[1]) > 1:
        writing_results(sys.argv[1])
    else:
        print("Error : you must specify a folder were the files will be created when you call this script !")
        print("Example : ")
        print("%(prog)s myfolder")
        exit(1)
