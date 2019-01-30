##########################
#   Description
##########################
#
# This program allow to retrieve every codon, anticodon, amino-acid from a set of exon given by the user
# Exons are designated by a name (give by the user), their chromosome and their chromosomal coordinates
# It also give a file containing the original input, the genomic and cds sequence of the exon

##########################
#   Imports
##########################

import mysql.connector  # allows the program to connect to the fasterDB database
import argparse  # for the parser
from exon_class import *  # allows the creation of an exon object
import commands
from dictionnary import dic_sequence  # useful to retrieve all the CDS of each exon in fasterDB
from multiprocessing import Pool  # allows to split the generation of control set in multiple processes
from exon_set_class import *  # allows the construction of a control exon
import random  # allows the random selection of control exons
from list_exon import *  # useful to manage a list of exons
import copy  # to allow deep copy of the table hsapiens_exonsstatus_improved
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector
from graphic import *
from composition_summary import *
from enrichiment_report_maker import writing_enrichment_report_file
from enrichiment_report_maker import writing_query_results_file
from Bio import SeqIO
import config

##########################
# Functions
##########################


def connection():
    """
    :return: an object that contains all the information you need to connect to fasterDB
    """
    cnx = mysql.connector.connect(user=config.user, password=config.password, host=config.host, database=config.database,
                                  buffered=True)
    return cnx


def retrieve_matching_exon(cnx, input_chromosome_number, input_chromosomal_coordinates_tuple, input_name):
    """
    :param cnx: the information necessary to connect to fasterDB
    :param input_chromosome_number: the chromosome of the exon of interest
    :param input_chromosomal_coordinates_tuple: a tuple containing the chromosomal coordinates of the exon
    :param input_name: the name given to the exon in the input
    :return: the list of all exons that are matching with those specifics chromosomal coordinates
    """

    cursor = cnx.cursor()
    query = ("""
            SELECT DISTINCT  id_gene, pos_sur_gene, start_sur_chromosome, end_sur_chromosome,
            NEW_cds_start_on_chromosome, NEW_cds_end_on_chromosome,exon_types, strand
            FROM  hsapiens_exonsstatus_improved
            WHERE ((end_sur_chromosome >= """ + str(input_chromosomal_coordinates_tuple[0]) + """
            AND end_sur_chromosome <= """ + str(input_chromosomal_coordinates_tuple[1]) + """)
            OR (start_sur_chromosome <= """ + str(input_chromosomal_coordinates_tuple[1]) + """
            AND end_sur_chromosome >= """ + str(input_chromosomal_coordinates_tuple[1]) + """))
            AND chromosome = \"""" + str(input_chromosome_number) + """\" ; """)
    cursor.execute(query)
    result = ListExon()
    for exon in cursor:
        result.exon_list.append(ExonClass([exon, input_name, input_chromosome_number]))
    if len(result.exon_list) == 0:  # no exon retrieved
        return None
    if len(result.exon_list) > 1:
        for i in range(len(result.exon_list)):
            result.exon_list[i].find_matching_status(len(result.exon_list))
    elif len(result.exon_list) == 1:
        result.exon_list[0].find_matching_status(len(result.exon_list))
    return result


def stop_message(message):
    """
    :param message: a message to display to the user
    If the sanitization check (see the function below) detects a problem in the input file, this function is called
    to inform the user that the program will stop
    """
    print "Input did not pass sanitization check."
    print message
    print "Terminating..."
    exit(1)


def sanitization_check(input_file):
    """
    :param input_file: the input file containing the exons that the user wants to study
    If an input line is not ok, nothing is done, the program stops
    :return: the number of line of the request
    """
    request = open(input_file, "r")
    line = request.readline()
    line = line.replace(" ", "")
    counter = 1
    while line:
        line = line.split("\t")

        if len(line) != 4:
            stop_message("Invalid number of columns in line " + str(counter))
        if str(line[1]) not in ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13",
                                "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y"]:
            stop_message("The chromosome number in line " + str(counter) + " is invalid")
        try:
            int(line[2])
            int(line[3])
        except ValueError:
            stop_message("Invalid chromosomal coordinates in line " + str(counter))
            exit()
        line = request.readline()
        counter += 1
    request.close()
    return counter


def sanitization_check_fasta(input_file):
    """
    :param input_file: the fasta file containing the sequence that the user wants to study
    If an input line is not ok, nothing is done, the program stops
    :return: the number of line of the request
    """
    i = 1
    for record in SeqIO.parse(input_file, "fasta"):
        for letter in record.seq:
            if letter not in "ATCGatgc":
                stop_message("Wrong sequence at sequence " + str(i))
        i += 1
    return i


def checking_output_path(output_folder):
    """
    :param output_folder: the path of the output folder
    :return: the final path that will be used (in case that the user has specified a wrong directory or is the user
    didn't specified any output directory
    """
    if output_folder[-1] != "/":
        output_folder += "/"
    if output_folder == "./" or os.path.isdir(output_folder) is False:
        if not os.path.isdir(output_folder):
            print "The folder " + output_folder + " doesn't exist ! "
            print "the output will be placed in the folder 'Exon_analysis' of your current directory !"
        if not os.path.isdir("./Exon_analysis"):
            os.mkdir("Exon_analysis")
        out_path = "./Exon_analysis/"
    else:
        out_path = output_folder
    return out_path


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


def checking_procedure(full_exon_list, complete_exon, correspondence, duplication):
    """
    :param full_exon_list: A list of the input exon
    :param complete_exon: A boolean value indicating if the user only wants complete exon
    :param correspondence: a boolean value indicating if the user wants the gene name and the exon position on this gene
    :param duplication: (boolean) True if we want to accept the duplication of a gene interval, False if we want to
     remove a duplicated interval
    corresponds to the user's exon name. The user exon name should in that case correspond to a concatenation of the
    gene name followed by an underscore and the exon position on this gene
    """

    if complete_exon is True and correspondence is False:
        i = 0
        while i < len(full_exon_list.exon_list):
            if not isinstance(full_exon_list.exon_list[i], list):
                if full_exon_list.exon_list[i].matching[1] != 100.0 or full_exon_list.exon_list[i].matching[2] != 100.0:
                    full_exon_list.exon_list[i] = [full_exon_list.exon_list[i].exon_name,
                                                   "Input_not_pass_user's_filters"]
                    i -= 1
            i += 1

    elif complete_exon is True and correspondence is True:
        i = 0
        while i < len(full_exon_list.exon_list):
            if not isinstance(full_exon_list.exon_list[i], list):
                if full_exon_list.exon_list[i].matching[1] != 100.0 or full_exon_list.exon_list[i].matching[2] != 100.0\
                        or full_exon_list.exon_list[i].exon_name.split("_")[0] != full_exon_list.exon_list[i].gene_name\
                        or full_exon_list.exon_list[i].exon_name.split("_")[1] != full_exon_list.exon_list[
                            i].exon_number:
                    full_exon_list.exon_list[i] = [full_exon_list.exon_list[i].exon_name,
                                                   "Input_not_pass_user's_filters"]
                    i -= 1
            i += 1
    elif complete_exon is False and correspondence is True:
        i = 0
        while i < len(full_exon_list.exon_list):
            if not isinstance(full_exon_list.exon_list[i], list):
                if full_exon_list.exon_list[i].exon_name.split("_")[0] != full_exon_list.exon_list[i].gene_name \
                        or full_exon_list.exon_list[i].exon_name.split("_")[1] != full_exon_list.exon_list[
                            i].exon_number:
                    full_exon_list.exon_list[i] = [full_exon_list.exon_list[i].exon_name,
                                                   "Input_not_pass_user's_filters"]
                    i -= 1
            i += 1

    if not duplication:
        exon_intervals = []
        for i in range(len(full_exon_list.exon_list)):
            if not isinstance(full_exon_list.exon_list[i], list):
                cur_interval = str(full_exon_list.exon_list[i].chr) + ":" + \
                               str(full_exon_list.exon_list[i].exon_start) + "-" + \
                               str(full_exon_list.exon_list[i].exon_end) + ";" + \
                               str(full_exon_list.exon_list[i].strand)
                if cur_interval in exon_intervals:
                    full_exon_list.exon_list[i] = [full_exon_list.exon_list[i].exon_name,
                                                   "Input_not_pass_user's_filters"]
                else:
                    exon_intervals.append(cur_interval)


def execute_request_fasta(input_file, out_path):
    """
    :param input_file: the file containing sequence we want to study
    :param out_path: the folder were the output files will be dropped
    :return: An exon_list objects containing exon objects : i.e. the exon given by the user
    """
    input_size = sanitization_check_fasta(input_file)  # checks if the input is correct
    input_content = [["fasta_name", "sequence"]]
    full_exon_list = ListExon()
    counter = 0
    for record in SeqIO.parse(input_file, "fasta"): # reading the input file
        if input_size > 300:
            print_progress(counter, input_size, "Processing input")
        counter += 1
        name_exon = str(record.id.split("|")[0])
        exon = ExonClass([name_exon, str(record.seq)], is_fasta=True)
        if len(str(record.seq)) <= 20:
            seq =  str(record.seq)
        else:
            seq = str(record.seq)[0:10] + "..." + str(record.seq)[-10:-1] + str(record.seq)[-1]
        input_content.append([name_exon, seq])
        exon.found_codon_anticodon_amino_acid_and_aa_nature()
        full_exon_list.exon_list.append(exon)

    # writing output files
    writing_query_results_file(out_path, input_content, full_exon_list)
    return full_exon_list, input_content


def execute_request_exons(input_file, out_path, complete_exon, correspondence, duplication):
    """
    :param input_file: the file containing the input of the exon
    :param out_path: the folder were the output files will be dropped
    :param complete_exon: true if the user only wants complete exon
    Run the programs : calculate for each interval (that normally corresponds to an exon) all exons
    that are in this interval, their sequences (genomic, cds peptide) and their composition (in codons,
    corresponding anticodons and amino acids)
    :param correspondence: a boolean value indicating if the user wants the gene name and the exon position on this gene
    corresponds to the user's exon name. The user exon name should in that case correspond to a concatenation of the
    gene name followed by an underscore and the exon position on this gene
    :param duplication: (boolean) True if we want to accept the duplication of a gene interval, False if we want to
     remove a duplicated interval
    :return: An exon_list objects containing exon objects : i.e. the exon given by the user
    """
    input_size = sanitization_check(input_file)  # checks if the input is correct
    request = open(input_file, "r")
    input_content = [["Exon_name", "Chromosome", "Start", "Stop"]]
    line = request.readline()
    cnx = connection()  # loading all the information to connect to fasterDB
    full_exon_list = ListExon()
    counter = 0
    while line:  # reading the input file
        if input_size > 300:
            print_progress(counter, input_size, "Processing input")
            counter += 1
        line = line.replace("\n", "").split("\t")
        input_content.append(line)  # copying each line a list for the query_results file
        line[2] = int(line[2])
        line[3] = int(line[3])
        list_of_exon = retrieve_matching_exon(cnx, line[1], line[2:4], line[0])
        if list_of_exon is not None:  # if the interval given by the user match exons we :
            for i in range(len(list_of_exon.exon_list)):
                list_of_exon.exon_list[i].retrieve_gene_name(cnx)  # retrieve the name of the genes that contains it
                list_of_exon.exon_list[i].retrieve_exon_sequences(cnx)  # retrieve its sequences
                list_of_exon.exon_list[i].calculate_exon_coverage_on_input_and_input_coverage_on_exon(line[2:4])
                list_of_exon.exon_list[i].resize_cds_sequence(
                    line[2:4])  # resize the sequences if the user only select a part of the exon
                # compute all the exon, corresponding anticodons, amino acids within the exon
                list_of_exon.exon_list[i].found_codon_anticodon_amino_acid_and_aa_nature()
            # we increment those new exon with no lacking information in a list containing all exons
            full_exon_list.exon_list += list_of_exon.exon_list
        else:
            list_of_exon = [line[0], "request did not match any fasterDB coding exon"]
            full_exon_list.exon_list.append(
                list_of_exon)  # we increment those new exon with no lacking information in a list containing all exons

        line = request.readline()

    checking_procedure(full_exon_list, complete_exon, correspondence, duplication)

    request.close()

    # writing output files
    writing_query_results_file(out_path, input_content, full_exon_list)
    return full_exon_list, input_content


def execute_request(input_file, out_path, complete_exon, correspondence, duplication):
    if input_file.split(".")[-1] in ["fasta", "fa"]:
        full_exon_list, input_content = execute_request_fasta(input_file, out_path)
    else:
        full_exon_list, input_content = execute_request_exons(input_file, out_path, complete_exon, correspondence, duplication)
    return full_exon_list, input_content


def concatenate(exon_list1, exon_list2):
    """
    Function to perform concatenation of the up exon list and the down exon list given by the user
    :param exon_list1: (list of exon instance), a list of exon
    :param exon_list2: (list of exon instance), another list of exon
    :return: (list of exon instance) the concatenation of the  exon_list1 and exon_list2
    """
    new_list = ListExon()
    for i in range(len(exon_list1.exon_list)):
        new_list.exon_list.append(exon_list1.exon_list[i])
    for i in range(len(exon_list2.exon_list)):
        new_list.exon_list.append(exon_list2.exon_list[i])
    return new_list


def concatenate_up_and_down_exon_list(output_folder, up_exon_list, down_exon_list, input_content):
    """
    Function to write some basic information on the combination of the down and the up exon list given by the user
    :param output_folder: (string) , the folder where the files (about the up and down exon list) will be stored
    :param up_exon_list: (list of exon instance) the user list of the up-regulated exon
    :param down_exon_list:(list of exon instance) the user list of the down-regulated exon
    :param input_content: (list of list of list) the concatenate content of the input files (exon up and down) given
    by the user. Each sublist corresponds to a line in the input files given by the user
    :return:(list of exon instance) the combination of the up and down exon list
    """
    full_exon_list = concatenate(up_exon_list, down_exon_list)
    writing_query_results_file(output_folder, input_content, full_exon_list)
    return full_exon_list


def calculate_interval(exon_size):
    """
    :param exon_size: the size of an exon
    :return: an interval near to the exon size
    This interval of value is useful when the user wants the control exons to have sizes close to his set of exons
    """
    if exon_size < 50:
        return max(exon_size / 3, 3), max(exon_size * 3, 50)
    elif exon_size > 350:
        return 190, 50000
    else:
        return exon_size / 2, exon_size * 2


def set_creator(set_size, user_exon_size, size_control, dic_sequences, query, penalty_size):
    """
    Function that performs the creation of a control set of exon
    :param set_size: (int) the size of the exon set : i.e. the number of exon in the set
    :param user_exon_size: (list of int) a list containing the sizes of the exons of the user
    :param size_control: (boolean) a boolean value indicating if the user wants (True) or not (False) that the exons of
    the control sets share a similar length than his exons
    :param dic_sequences: (dictionary of string) a dictionary containing all the cds of exons from a type
    (ACE, LCE, FCE, CCE, ALL) in fasterDB database
    :param query: (list of tuple) a tuple containing all the exon from a type (ACE, LCE, FCE, CCE, ALL) in the fasterDB
    database
    :param penalty_size: (int) size in nucleotide of exon that will be considered as short when frequencies of codon,
    amino acid... are calculated and thus will have less weight than the others
    :return: three dictionary : one containing the frequency of each codon in the control set
    one containing the frequency of each amino acid in the control set
    one containing the frequency of different group of amino acid based on their side chain properties (schain group)
    one containing the frequency of different group of amino acid based on their hydrophilic/phobic properties
    (hydro group)
    one containing the frequency of different group of amino acid based on their polarity properties (polar group)
    one containing the frequency of different group of amino acid based on their charge properties (charge group)
    one containing the frequency of different group of amino acid based on miscellaneous properties (misc group)
    one containing the frequency of different group of amino acid based on their chemical properties (chemical group)
    one containing the frequency of different group of amino acid based on their structural properties(structural group)
    """
    query_result = copy.deepcopy(query)
    set_result = ListExon()
    for i in range(set_size):
        value = random.randint(1, len(query_result) - 1)
        if size_control is False:
            while int(query_result[value][2]) < 3:
                value = random.randint(1, len(query_result) - 1)
        else:
            interval = calculate_interval(user_exon_size[i])
            while int(query_result[value][2]) < interval[0] or int(query_result[value][2]) > interval[1]:
                value = random.randint(1, len(query_result) - 1)
        set_result.exon_list.append(ExonSetClass(query_result[value]))
        del (query_result[value])
    for i in range(len(set_result.exon_list)):
        set_result.exon_list[i].retrieve_exon_sequences(dic_sequences)
        set_result.exon_list[i].found_codon_anticodon_amino_acid_and_aa_nature()
    codon_frequency = set_result.weighted_frequency_calculator(penalty_size / 3)
    aa_frequency = set_result.amino_acid_frequency_calculator(penalty_size / 3)
    ft_frequency = set_result.weighted_feature_calculator(penalty_size / 3)
    ratio_frequency = set_result.feature_ratio(penalty_size / 3)
    opposed_ratio_frequency = set_result.opposed_feature_ratio(penalty_size / 3)
    nucleic_acid_frequency = set_result.nucleic_acid_calculator()
    nt_pos_frequency = set_result.combine_nt_pos_frequencies()
    dnt_frequency = set_result.dinucleotide_calculator()
    hexa_frequency = set_result.hexanucleotide_calculator(penalty_size)
    diaa_frequency = set_result.diaa_calculator(penalty_size / 3)
    return codon_frequency, aa_frequency, ft_frequency, ratio_frequency, opposed_ratio_frequency, \
           nucleic_acid_frequency, nt_pos_frequency, dnt_frequency, hexa_frequency, diaa_frequency


def sets_creator(set_size, set_number, user_exon_size, size_control, dic_sequences, query_result, progression,
                 penalty_size):
    """
    Function to perform the creation of all the control sets of exon
    :param set_size: (int) the size of the exon set. i.e. the number of exon in the set
    :param set_number: (int) the number of sets to create
    :param user_exon_size: (list of int) a list containing the sizes of the exons of the user
    :param size_control: (boolean) a boolean value indicating if the user wants (True) or not (False) that the exons of
    the control sets share a similar length than his exons
    :param dic_sequences: (dictionary of string) a dictionary containing all the cds of exons from a type
    (ACE, LCE, FCE, CCE, ALL) in fasterDB database
    :param query_result: (list of tuple) a tuple containing all the exon from a type (ACE, LCE, FCE, CCE, ALL) in the
    fasterDB database
    :param progression: (boolean) true if you and to see the progression
    :param penalty_size: (int) size in nucleotide of exon that will be considered as short when frequencies of codon,
    amino acid... are calculated and thus will have less weight than the others
    :return: a dictionary containing the codon, amino acid, and some other group
    (schain, hydro, charge, polarity, misc, chemical, structural) frequencies of the set_number sets.
    Thus for each codon/amino_acid/nature set_number frequencies are available, one for each set
    """
    codon_frequencies = dict()
    aa_frequencies = dict()
    ft_frequencies = dict()
    ratio_frequencies = dict()
    opposed_ratio_frequencies = dict()
    nucleic_acid_frequencies = dict()
    nt_pos_frequencies = dict()
    dnt_frequencies = dict()
    hexa_frequencies = dict()
    diaa_frequencies = dict()
    codon_frequency, aa_frequency, ft_frequency, ratio_frequency, opposed_ratio_frequency, \
    nucleic_acid_frequency, nt_pos_frequency, dnt_frequency, hexa_frequency, diaa_frequency = \
        set_creator(set_size, user_exon_size, size_control, dic_sequences, query_result, penalty_size)
    for codon in codon_frequency.keys():
        codon_frequencies[codon] = list()
        codon_frequencies[codon].append(codon_frequency[codon])
    for aa in aa_frequency.keys():
        aa_frequencies[aa] = list()
        aa_frequencies[aa].append(aa_frequency[aa])
    for ft in ft_frequency.keys():
        ft_frequencies[ft] = list()
        ft_frequencies[ft].append(ft_frequency[ft])
    for ftr in ratio_frequency.keys():
        ratio_frequencies[ftr] = list()
        ratio_frequencies[ftr].append(ratio_frequency[ftr])
    for ftor in opposed_ratio_frequency.keys():
        opposed_ratio_frequencies[ftor] = list()
        opposed_ratio_frequencies[ftor].append(opposed_ratio_frequency[ftor])
    for nt in nucleic_acid_frequency.keys():
        nucleic_acid_frequencies[nt] = list()
        nucleic_acid_frequencies[nt].append(nucleic_acid_frequency[nt])
    for ntp in nt_pos_frequency.keys():
        nt_pos_frequencies[ntp] = list()
        nt_pos_frequencies[ntp].append(nt_pos_frequency[ntp])
    for dnt in dnt_frequency.keys():
        dnt_frequencies[dnt] = list()
        dnt_frequencies[dnt].append(dnt_frequency[dnt])
    for hexa in hexa_frequency.keys():
        hexa_frequencies[hexa] = list()
        hexa_frequencies[hexa].append(hexa_frequency[hexa])
    for diaa in diaa_frequency.keys():
        diaa_frequencies[diaa] = list()
        diaa_frequencies[diaa].append(diaa_frequency[diaa])
    for i in range(set_number - 1):
        if progression is True:
            print_progress(i + 2, set_number, "Creating control sets")
        codon_frequency, aa_frequency, ft_frequency, ratio_frequency, opposed_ratio_frequency, \
        nucleic_acid_frequency, nt_pos_frequency, dnt_frequency, hexa_frequency, diaa_frequency = \
            set_creator(set_size, user_exon_size, size_control, dic_sequences, query_result, penalty_size)
        for codon in codon_frequency.keys():
            codon_frequencies[codon].append(codon_frequency[codon])
        for aa in aa_frequency.keys():
            aa_frequencies[aa].append(aa_frequency[aa])
        for ft in ft_frequency.keys():
            ft_frequencies[ft].append(ft_frequency[ft])
        for ftr in ratio_frequency.keys():
            ratio_frequencies[ftr].append(ratio_frequency[ftr])
        for ftor in opposed_ratio_frequency.keys():
            opposed_ratio_frequencies[ftor].append(opposed_ratio_frequency[ftor])
        for nt in nucleic_acid_frequency.keys():
            nucleic_acid_frequencies[nt].append(nucleic_acid_frequency[nt])
        for ntp in nt_pos_frequency.keys():
            nt_pos_frequencies[ntp].append(nt_pos_frequency[ntp])
        for dnt in dnt_frequency.keys():
            dnt_frequencies[dnt].append(dnt_frequency[dnt])
        for hexa in hexa_frequency.keys():
            hexa_frequencies[hexa].append(hexa_frequency[hexa])
        for diaa in diaa_frequency.keys():
            diaa_frequencies[diaa].append(diaa_frequency[diaa])
    return codon_frequencies, aa_frequencies, ft_frequencies, ratio_frequencies, opposed_ratio_frequencies, \
           nucleic_acid_frequencies, nt_pos_frequencies, dnt_frequencies, hexa_frequencies, diaa_frequencies


def calculate_index(value, list_of_values):
    """
    :param value: a value corresponding to a codon frequency
    :param list_of_values: a list of codon frequencies
    :return: an empirical "p-value" calculated by the position of the value "value" in the list_of_value
    sorted in the ascending order
    """
    list_of_values.sort()
    val1 = -1
    val2 = -1
    for i in range(len(list_of_values)):  # proportion of values in the list_of_values bigger than the value "value"
        if value <= list_of_values[i]:
            val1 = float(len(list_of_values) - i) / len(list_of_values)
            break  # as soon as we find the index, we quit the loop
    for i in range(len(list_of_values)):  # proportion of values in the list_of_values smaller than the value "value"
        if value < list_of_values[i]:
            val2 = float(i) / len(list_of_values)
            break  # as soon as we find the index, we quit the loop
    if val2 == -1 and list_of_values[-1] == value:
        return 1
    if val1 != -1 and val2 != -1:
        return min(min(val1, val2), 1)
    return 0


def calculate_enrichment(control_frequencies, interest_frequencies):
    """
    :param control_frequencies: frequencies of the control sets
    :param interest_frequencies: frequency of each codon in the user set of exons
    :return: a dictionary containing the "p-value" of each codon
    """
    dic_p_val = {}
    for codon_or_amino_acid_or_nature in control_frequencies.keys():
        dic_p_val[codon_or_amino_acid_or_nature] = calculate_index(interest_frequencies[codon_or_amino_acid_or_nature],
                                                                   control_frequencies[codon_or_amino_acid_or_nature])
    return dic_p_val


def sql_request(cnx, exon_type):
    """
    :param cnx: (mysql.connector instance object) all the information you need to connect to fasterDB
    :param exon_type: (string) the type of the exon (ACE, FCE, LCE, CCE, ALL)
    :return: a list of tuple corresponding to all the exon in fasterDB and their chromosomal coordinates
    """
    cursor = cnx.cursor()
    if exon_type != "ALL":
        query = ("""
                    SELECT DISTINCT  id_gene, pos_sur_gene, NEW_cds_end_on_chromosome-NEW_cds_start_on_chromosome+1
                    FROM  hsapiens_exonsstatus_improved
                    WHERE exon_types LIKE \"%""" + str(exon_type) + """%\" ;""")
    else:
        query = ("""
                    SELECT DISTINCT  id_gene, pos_sur_gene, NEW_cds_end_on_chromosome-NEW_cds_start_on_chromosome+1
                    FROM  hsapiens_exonsstatus_improved""")
    cursor.execute(query)
    query_result = cursor.fetchall()
    return query_result


def combine_a_result(list_of_results_from_the_threads, number_of_result):
    """
    :param number_of_result: (int) the number of the frequency dictionary to combine
    :param list_of_results_from_the_threads: a list of dictionaries, 3 dictionary per threads
    :return: the combined result of the threads for one type of frequency (codon, amino_acid, amino acid nature,
    amino acid importance, amino acid metabolism)
    """
    control_frequencies = dict()
    for key in list_of_results_from_the_threads[0][number_of_result].keys():
        frequency_completed = list()
        for i in range(len(list_of_results_from_the_threads)):
            frequency_completed += list_of_results_from_the_threads[i][number_of_result][key]
        control_frequencies[key] = frequency_completed
    return control_frequencies


def combine_all_results(list_of_results_from_the_threads):
    """
    :param list_of_results_from_the_threads: a list of dictionaries, 3 dictionary per threads
    :return: the combined results of the threads
    """
    result = list()
    for i in range(len(list_of_results_from_the_threads[0])):
        result.append(combine_a_result(list_of_results_from_the_threads, i))
    return result


def multiprocesses(thread, set_size, set_number, size_of_the_user_exon, size_control, dic_sequences, query,
                   penalty_size):
    """
    :param thread: (int) the number of processor the user wants to use to compute the control sets of exons
    :param set_size: (int) the size of each control set. It is the same as the user set
    :param set_number: (int) the number of control sets the user wants to create
    :param size_of_the_user_exon: (list of int) The size (in nucleotides) for each exons in the user set
    :param size_control: (boolean) a boolean value indicating if the user wants (True) or not (False) that the exons of
    the control sets share a similar length than his exons
    :param dic_sequences: (dictionary of string) a dictionary containing all the cds of exons from a type
    (ACE, LCE, FCE, CCE, ALL) in fasterDB database
    :param query: (list of tuple) a tuple containing all the exon from a type (ACE, LCE, FCE, CCE, ALL) in the fasterDB
    database
    :param penalty_size: (int) size in nucleotide of exon that will be considered as short when frequencies of codon,
    amino acid... are calculated and thus will have less weight than the others
    :return: the combine results (frequency of each codon, amino acid and their nature ofr each control sets) compute
    by the processor(s)
    """
    pool = Pool()
    process = list()
    for i in range(thread):
        if i == 0:
            process.append(pool.apply_async(sets_creator, [set_size, set_number / thread, size_of_the_user_exon,
                                                           size_control, dic_sequences, query, True, penalty_size]))
        else:
            process.append(pool.apply_async(sets_creator,
                                            [set_size, set_number / thread, size_of_the_user_exon,
                                             size_control, dic_sequences, query, False, penalty_size]))
    my_res = list()
    for i in range(len(process)):
        my_res.append(process[i].get(timeout=None))
    return combine_all_results(my_res)


def write_values_of_all_control_sets(control_dic, name_dic, output):
    """
    :param control_dic: (dictionary of list of float)
    :param name_dic: (string) the name of the dictionary
    :param output: (string ) the path where the result will be created
    """
    with open(output + name_dic + "_control_value.tsv", "w") as out_file:
        header = ""
        list_key = control_dic.keys()
        for key in list_key:
            header += key + "\t"
        header = header[:-1] + "\n"
        out_file.write(header)
        for i in range(len(control_dic[list_key[0]])):
            line = ""
            for key in list_key:
                line += str(control_dic[key][i]) + "\t"
            line = line[:-1] + "\n"
            out_file.write(line)


def main_program(full_exon_list, output_folder, enrichment, exon_type, set_number, size_control, thread, penalty_size,
                 name):
    """
    :param full_exon_list: (list of exon instance) an exon_list object containing exon object
    :param output_folder: (string) the output folder were the user want to drop off his result
    :param enrichment: (boolean) a boolean value indicating if the user want to do an enrichment analysis
    :param exon_type: (string) the type of the exon of the control sets
    :param set_number: (int)  the number of control sets to create
    :param size_control: (boolean) a boolean value indicating if the user wants (True) or not (False) that the exons of
    the control sets share a similar length than his exons
    :param thread: (int) the number of processor the user wants to use to compute the control sets of exons
    :param penalty_size: (int) size in nucleotide of exon that will be considered as short when frequencies of codon,
    amino acid... are calculated and thus will have less weight than the others
    :param name: (string) corresponding to the name of the set of exon that will figure in the graphics produced by the
    program
    """
    cnx = connection()
    set_size = len(full_exon_list.exon_list)
    if enrichment is False:
        del full_exon_list
        exit()
    else:

        #########################################
        # Test frequencies  ########
        freq_codon_per_exon = full_exon_list.weighted_exon_frequency_calculator(round(penalty_size / 3))
        freq_aa_per_exon = full_exon_list.amino_acid_exon_frequency_calculator(round(penalty_size / 3))

        # writing detail of codon and amino acid
        out_path = checking_output_path(output_folder)
        file_codon = open(out_path + "frequency_codon_foreach_exon.info", "w")
        for key in freq_codon_per_exon.keys():
            file_codon.write(key + "\t" + str(freq_codon_per_exon[key]).replace("[", "").replace("]", "") + "\n")
        file_codon.close()
        file_aa = open(out_path + "frequency_aa_foreach_exon.info", "w")
        for key in freq_aa_per_exon.keys():
            file_aa.write(key + "\t" + str(freq_aa_per_exon[key]).replace("[", "").replace("]", "") + "\n")
        file_aa.close()

        ##########################################

        dic_sequences = dic_sequence(cnx, exon_type)
        query = sql_request(cnx, exon_type)
        size_of_the_user_exon = full_exon_list.return_the_size_of_each_exon_in_the_list()

        nucleotides_exons_frequencies = full_exon_list.mean_exon_nucleotide_calculator()
        nucleotide_frequencies = full_exon_list.weighted_nucleotide_calculator(round(penalty_size))

        codon_frequencies = full_exon_list.weighted_frequency_calculator(round(penalty_size / 3))
        codon_frequencies_5p = full_exon_list.weighted_frequency_calculator(round(penalty_size / 3), extr = "5p")
        codon_frequencies_3p = full_exon_list.weighted_frequency_calculator(round(penalty_size / 3), extr = "3p")

        aa_frequencies = full_exon_list.amino_acid_frequency_calculator(round(penalty_size / 3))
        aa_frequencies_5p = full_exon_list.amino_acid_frequency_calculator(round(penalty_size / 3), extr = "5p")
        aa_frequencies_3p = full_exon_list.amino_acid_frequency_calculator(round(penalty_size / 3), extr = "3p")

        ft_frequency = full_exon_list.weighted_feature_calculator(penalty_size / 3,)
        ft_frequency_5p = full_exon_list.weighted_feature_calculator(penalty_size / 3, extr = "5p")
        ft_frequency_3p = full_exon_list.weighted_feature_calculator(penalty_size / 3, extr = "3p")

        ftr_frequency = full_exon_list.feature_ratio(penalty_size / 3,)
        ftr_frequency_5p = full_exon_list.feature_ratio(penalty_size / 3, extr="5p")
        ftr_frequency_3p = full_exon_list.feature_ratio(penalty_size / 3, extr="3p")

        ftor_frequency = full_exon_list.opposed_feature_ratio(penalty_size / 3)
        ftor_frequency_5p = full_exon_list.opposed_feature_ratio(penalty_size / 3, extr="5p")
        ftor_frequency_3p = full_exon_list.opposed_feature_ratio(penalty_size / 3, extr="3p")

        nucleic_acid_frequency = full_exon_list.nucleic_acid_calculator()
        nucleic_acid_frequency_5p = full_exon_list.nucleic_acid_calculator(extr = "5p")
        nucleic_acid_frequency_3p = full_exon_list.nucleic_acid_calculator(extr = "3p")

        ntp_frequency = full_exon_list.combine_nt_pos_frequencies()
        ntp_frequency_5p = full_exon_list.combine_nt_pos_frequencies(extr="5p")
        ntp_frequency_3p = full_exon_list.combine_nt_pos_frequencies(extr="3p")

        dnt_frequency = full_exon_list.dinucleotide_calculator()
        dnt_frequency_5p = full_exon_list.dinucleotide_calculator(extr = "5p")
        dnt_frequency_3p = full_exon_list.dinucleotide_calculator(extr = "3p")

        hexa_frequency = full_exon_list.hexanucleotide_calculator(penalty_size)
        hexa_frequency_5p = full_exon_list.hexanucleotide_calculator(penalty_size, extr = "5p")
        hexa_frequency_3p = full_exon_list.hexanucleotide_calculator(penalty_size, extr = "3p")

        diaa_frequency = full_exon_list.diaa_calculator(penalty_size / 3)
        diaa_frequency_5p = full_exon_list.diaa_calculator(penalty_size / 3 , extr = "5p")
        diaa_frequency_3p = full_exon_list.diaa_calculator(penalty_size / 3 , extr = "3p")

        if thread != 1:
            control_frequencies_codon, control_frequencies_aa, control_ft_frequencies, \
            control_ftr_frequencies, control_ftor_frequencies, control_nucleic_acid_frequencies, \
            control_ntp_frequencies, control_dnt_frequencies, control_hexa_frequencies, control_diaa_frequencies = \
                multiprocesses(thread, set_size, set_number, size_of_the_user_exon, size_control, dic_sequences, query,
                           penalty_size)
        else:
            control_frequencies_codon, control_frequencies_aa, control_ft_frequencies, \
            control_ftr_frequencies, control_ftor_frequencies, control_nucleic_acid_frequencies, \
            control_ntp_frequencies, control_dnt_frequencies, control_hexa_frequencies, control_diaa_frequencies = \
                sets_creator(set_size, set_number, size_of_the_user_exon, size_control, dic_sequences, query, True,
                         penalty_size)

        dic_p_val_codon = calculate_enrichment(control_frequencies_codon, codon_frequencies)
        dic_p_val_aa = calculate_enrichment(control_frequencies_aa, aa_frequencies)
        dic_p_val_ft = calculate_enrichment(control_ft_frequencies, ft_frequency)
        dic_p_val_ftr = calculate_enrichment(control_ftr_frequencies, ftr_frequency)
        dic_p_val_ftor = calculate_enrichment(control_ftor_frequencies, ftor_frequency)
        dic_p_val_nt = calculate_enrichment(control_nucleic_acid_frequencies, nucleic_acid_frequency)
        dic_p_val_ntp = calculate_enrichment(control_ntp_frequencies, ntp_frequency)
        dic_p_val_dnt = calculate_enrichment(control_dnt_frequencies, dnt_frequency)
        dic_p_val_hexa = calculate_enrichment(control_hexa_frequencies, hexa_frequency)
        dic_p_val_diaa = calculate_enrichment(control_diaa_frequencies, diaa_frequency)

        write_values_of_all_control_sets(control_ntp_frequencies, "nucleotide_position", output_folder)

        dic_padjust_codon, dic_padjust_aa = \
            writing_enrichment_report_file(control_frequencies_codon, codon_frequencies, codon_frequencies_5p, codon_frequencies_3p, dic_p_val_codon,
                                           control_frequencies_aa, aa_frequencies, aa_frequencies_5p, aa_frequencies_3p, dic_p_val_aa,
                                           control_ft_frequencies, ft_frequency, ft_frequency_5p, ft_frequency_3p, dic_p_val_ft,
                                           control_ftr_frequencies, ftr_frequency, ftr_frequency_5p, ftr_frequency_3p, dic_p_val_ftr,
                                           control_ftor_frequencies, ftor_frequency, ftor_frequency_5p, ftor_frequency_3p, dic_p_val_ftor,
                                           control_nucleic_acid_frequencies, nucleic_acid_frequency, nucleic_acid_frequency_5p, nucleic_acid_frequency_3p, dic_p_val_nt,
                                           control_ntp_frequencies, ntp_frequency, ntp_frequency_5p, ntp_frequency_3p, dic_p_val_ntp,
                                           control_dnt_frequencies, dnt_frequency, dnt_frequency_5p, dnt_frequency_3p, dic_p_val_dnt,
                                           control_hexa_frequencies, hexa_frequency, hexa_frequency_5p, hexa_frequency_3p, dic_p_val_hexa,
                                           control_diaa_frequencies, diaa_frequency, diaa_frequency_5p, diaa_frequency_3p, dic_p_val_diaa,
                                           output_folder, set_number)
        create_graphs(codon_frequencies, aa_frequencies, dic_p_val_codon, dic_p_val_aa, dic_padjust_codon,
                      dic_padjust_aa, exon_type, out_path, name)

        if "up" in name or "down" in name:
            return codon_frequencies, aa_frequencies, nucleotide_frequencies, nucleotides_exons_frequencies


def control_summary_set_creator(cnx, penalty_size, exon_type):
    """
    :param cnx: (mysql.connector object) an object that contains all the information you need to connect to fasterDB
    :param penalty_size: (int) size in nucleotide of exon that will be considered as short when frequencies of codon,
    amino acid... are calculated and thus will have less weight than the others
    :param exon_type: (string) the type of the exon (ACE, FCE, LCE, CCE, ALL)
    :return: 7 dictionaries
    metaexon method : all the sequence in a exon set are concatenated. After that the frequency calculation is done
    on the meta-sequence obtained
    weighted method : all the calculation on codon/amino acid are done for each exon sequence and are then pondered
    by the number of exons in the exon set or by number of exons in the exon set and the exon size if the exon is small.
    An exon is considered as small is its size in smaller than the penalty_size value.
    for example is we have 3 exons ( ATG-ACC-ATG  ATG_ATC AAA-TTT) and the penalty size is 9 / 3 = 3 codons
    the frequency of ATG within this set of 3 exon will be : (0.5 + 0.5 * 2/3 + 0) / (1 + 2/3 + 1)
    1 - meta_codon_frequency : for each codon, gives their frequencies in all the exons "exon_types" of fasterDB
    (with the metaexon method - see list_exon)
    2 - meta_aa_frequency : for each amino acids, gives their frequencies in all the exons "exon_types" of fasterDB
    (with the metaexon method - see list_exon)
    3 - meta_codon_composition : gives the frequency of codons that are rich or poor in A, T, G, C, AC, AG, AT, CG, CT,
    GT in all the exon "exon_type" of fasterDB (with the metaexon method - see list_exon)
    4 - meta_codon_last_nt : gives the frequency of codons that end by A, T, G, C, in all the exon "exon_type" of
    fasterDB (with the metaexon method - see list_exon)
    5 - weighted_codon_frequency : for each codon, gives their frequencies in all the exons "exon_types" of fasterDB
    (with the weighted method - see list_exon)
    6 - weighted_aa_frequency : for each amino acids, gives their frequencies in all the exons "exon_types" of fasterDB
    (with the weighted method - see list_exon)
    7 - weighted_composition : gives the frequency of codons that are rich or poor in A, T, G, C, AC, AG, AT, CG, CT,
    GT in all the exon "exon_type" of fasterDB (with the weighted method - see list_exon)
    8 - weighted_last_nt : gives the frequency of codons that end by A, T, G, C, in all the exon "exon_type" of
    fasterDB (with the weighted method - see list_exon)
    """
    query_result = sql_request(cnx, exon_type)
    dic_sequences = dic_sequence(cnx, exon_type)
    set_result = ListExon()
    for i in range(len(query_result)):
        set_result.exon_list.append(ExonSetClass(query_result[i]))
    for i in range(len(set_result.exon_list)):
        set_result.exon_list[i].retrieve_exon_sequences(dic_sequences)
        set_result.exon_list[i].found_codon_amino_acid_and_aa_nature()
    print "calculating the dictionaries for " + str(exon_type) + " exons..."
    # meta_exon method for the recap
    meta_codon_frequency = set_result.meta_exon_codon_calculator()
    meta_aa_frequency = set_result.meta_exon_aa_calculator()
    meta_codon_composition = set_result.meta_exon_composition_analyser()
    meta_codon_last_nt = set_result.meta_exon_last_nt_proportion_calculator()
    # weighted exon method
    weighted_codon_frequency = set_result.weighted_frequency_calculator_recap(penalty_size / 3)
    weighted_aa_frequency = set_result.amino_acid_frequency_calculator_recap(penalty_size / 3)
    weighted_composition = set_result.weighted_exon_composition_recap(penalty_size / 3)
    weighted_last_nt = set_result.weighted_last_nt_frequency_calculator_recap(penalty_size / 3)
    return meta_codon_frequency, meta_aa_frequency, meta_codon_composition, meta_codon_last_nt, \
           weighted_codon_frequency, weighted_aa_frequency, weighted_composition, weighted_last_nt


def creating_composition_summary(output_folder, up_exon_list, down_exon_list, up_and_down_list, penalty_size):
    """
    Function that perform the creation of the composition_summary.xlsx file
    :param output_folder: (string) the path where the user wants to create the composition_summary.xlsx file
    :param up_exon_list: (list of exon instance) the list of exon up-regulated in a particular condition given
    by the user
    :param down_exon_list:(list of exon instance) the list of exon down-regulated in a particular condition given
    by the user
    :param up_and_down_list:(list of exon instance) the list of exon up and down regulated in a particular condition
    given by the user
    :param penalty_size: (int) size in nucleotide of exon that will be considered as short when frequencies of codon,
    amino acid... are calculated and thus will have less weight than the others
    """
    cnx = connection()
    # getting information on the three list of exons studied
    dic_up = getting_summary_dictionaries(up_exon_list, penalty_size)
    dic_down = [None] * 8
    dic_up_down = [None] * 8
    if down_exon_list is not None:
        dic_down = getting_summary_dictionaries(down_exon_list, penalty_size)
    if up_and_down_list is not None:
        dic_up_down = getting_summary_dictionaries(up_and_down_list, penalty_size)

    # getting information on the 2 control list of exons (CCE and ACE control list of exon)
    dic_ace = control_summary_set_creator(cnx, penalty_size, "ACE")
    dic_cce = control_summary_set_creator(cnx, penalty_size, "CCE")

    # creating the content for the sheet
    content, merging_values, merging_content = \
        get_content_first_sheet(dic_up[0], dic_up[1], dic_down[0], dic_down[1], dic_up_down[0], dic_up_down[1],
                                dic_ace[0], dic_ace[1], dic_cce[0], dic_cce[1])
    content2, merging_values2, merging_content2 = \
        get_content_second_sheet(dic_up[2], dic_up[3], dic_down[2], dic_down[3], dic_up_down[2], dic_up_down[3],
                                 dic_ace[2], dic_ace[3], dic_cce[2], dic_cce[3])
    content3, merging_values3, merging_content3 =\
        get_content_first_sheet(dic_up[4], dic_up[5], dic_down[4], dic_down[5], dic_up_down[4], dic_up_down[5],
                                dic_ace[4], dic_ace[5], dic_cce[4], dic_cce[5])
    content4, merging_values4, merging_content4 = \
        get_content_second_sheet(dic_up[6], dic_up[7], dic_down[6], dic_down[7], dic_up_down[6], dic_up_down[7],
                                 dic_ace[6], dic_ace[7], dic_cce[6], dic_cce[7])

    write_xls(output_folder, content, merging_values, merging_content, content2, merging_values2, merging_content2,
              content3, merging_values3, merging_content3, content4, merging_values4, merging_content4)


def creating_recap_graphics(output, up_exon_list, down_exon_list, up_down_exon_list):
    """
    Function that creates the summary figures
    :param output: (string) the path where the figure folder will be created
    :param up_exon_list: (list of exon instance) the list of exon up-regulated in a particular condition given
    by the user
    :param down_exon_list: (list of exon instance) the list of exon down-regulated in a particular condition given
    by the user
    :param up_down_exon_list:(list of exon instance) the list of exon up and down regulated in a particular condition
    given by the user
    """

    # create a working directory
    os.mkdir(output + "recap_graphics/")
    # for the boxplot
    os.mkdir(output + "recap_graphics/boxplot/")
    # One for the codon graphics
    os.mkdir(output + "recap_graphics/boxplot/codon/")
    # One for the amino acid graphics
    os.mkdir(output + "recap_graphics/boxplot/amino_acid/")
    # for the barplot
    os.mkdir(output + "recap_graphics/barplot/")
    # One for the codon graphics
    os.mkdir(output + "recap_graphics/barplot/codon/")
    # One for the amino acid graphics
    os.mkdir(output + "recap_graphics/barplot/amino_acid/")

    # creating boxplot for codons
    print "Creating boxplot summary figure for codon..."

    for codon in codon2aminoAcid.keys():
        # obtaining the value we need to construct the graphics
        up_prop_list = up_exon_list.get_exon_proportion_of_codon(codon)
        down_prop_list = None
        if down_exon_list is not None:
            down_prop_list = down_exon_list.get_exon_proportion_of_codon(codon)
        up_down_prop_list = None
        if up_down_prop_list is not None:
            up_down_prop_list = up_down_exon_list.get_exon_proportion_of_codon(codon)

        # creating the graphics...
        boxplot_composition(output + "recap_graphics/boxplot/codon/", up_prop_list, down_prop_list, up_down_prop_list,
                            codon)

    print "Creating boxplot for amino_acid..."
    for aa in amino_acid2nature.keys():
        # obtaining the value we need to construct the graphics
        up_prop_list = up_exon_list.get_exon_proportion_of_amino_acid(aa)
        down_prop_list = None
        if down_exon_list is not None:
            down_prop_list = down_exon_list.get_exon_proportion_of_amino_acid(aa)
        up_down_prop_list = None
        if up_down_prop_list is not None:
            up_down_prop_list = up_down_exon_list.get_exon_proportion_of_amino_acid(aa)

        # creating the graphics...
        boxplot_composition(output + "recap_graphics/boxplot/amino_acid/", up_prop_list, down_prop_list,
                            up_down_prop_list, aa)
    # ------------------------- barplot
    print "Creating barplot summary figure for codon..."

    for codon in codon2aminoAcid.keys():
        # obtaining the value we need to construct the graphics
        up_prop_list = up_exon_list.get_exon_proportion_count_of_codon(codon)
        down_prop_list = None
        if down_exon_list is not None:
            down_prop_list = down_exon_list.get_exon_proportion_count_of_codon(codon)

        # creating the graphics...
        barplot_count_proportion(output + "recap_graphics/barplot/codon/", up_prop_list, down_prop_list, codon)
        fichier = open(output + "recap_graphics/barplot/up_count_codon" ,"a")
        fichier.write(codon + "\t" + str(up_prop_list[0]) + "\t" + str(up_prop_list[1]) + "\t" + str(up_prop_list[2]) + "\t" + str(up_prop_list[3]) + "\t" + str(up_prop_list[4]) + "\n")
        fichier.close()
        if down_exon_list is not None:
            fichier = open(output + "recap_graphics/barplot/down_count_codon" ,"a")
            fichier.write(codon + "\t" + str(down_prop_list[0]) + "\t" + str(down_prop_list[1]) + "\t" + str(down_prop_list[2]) + "\t" + str(down_prop_list[3]) + "\t" + str(down_prop_list[4]) + "\n")
            fichier.close()
    print "Creating boxplot for amino_acid..."
    for aa in amino_acid2nature.keys():
        # obtaining the value we need to construct the graphics
        up_prop_list = up_exon_list.get_exon_proportion_count_of_amino_acid(aa)
        down_prop_list = None
        if down_exon_list is not None:
            down_prop_list = down_exon_list.get_exon_proportion_count_of_amino_acid(aa)

        # creating the count graphics...
        barplot_count_proportion(output + "recap_graphics/barplot/amino_acid/", up_prop_list, down_prop_list, aa)
        fichier = open(output + "recap_graphics/barplot/up_count_aa" ,"a")
        fichier.write(aa + "\t" + str(up_prop_list[0]) + "\t" + str(up_prop_list[1]) + "\t" + str(up_prop_list[2]) + "\t" + str(up_prop_list[3]) + "\t" + str(up_prop_list[4]) + "\n")
        fichier.close()
        if down_exon_list is not None:
            fichier = open(output + "recap_graphics/barplot/down_count_aa" ,"a")
            fichier.write(aa + "\t" + str(down_prop_list[0]) + "\t" + str(down_prop_list[1]) + "\t" + str(down_prop_list[2]) + "\t" + str(down_prop_list[3]) + "\t" + str(down_prop_list[4]) +"\n")
            fichier.close()

        # creating the count graphics...
        barplot_count_proportion(output + "recap_graphics/barplot/amino_acid/", up_prop_list, down_prop_list, aa)
        fichier = open(output + "recap_graphics/barplot/up_prop_aa" ,"a")
        fichier.write(aa + "\t" + str(float(up_prop_list[0])/up_prop_list[4]) + "\t" + str(float(up_prop_list[1])/up_prop_list[4]) + "\t" + str(float(up_prop_list[2])/up_prop_list[4]) + "\t" + str(float(up_prop_list[3])/up_prop_list[4]) + "\n")
        fichier.close()
        if down_exon_list is not None:
            fichier = open(output + "recap_graphics/barplot/down_prop_aa" ,"a")
            fichier.write(aa + "\t" + str(float(down_prop_list[0])/down_prop_list[4]) + "\t" + str(float(down_prop_list[1])/down_prop_list[4]) + "\t" + str(float(down_prop_list[2])/down_prop_list[4]) + "\t" + str(float(down_prop_list[3])/down_prop_list[4]) +"\n")
            fichier.close()


def up_and_down_launcher(up_exon_list, down_exon_list, up_and_down_list, output_folder, enrichment, exon_type,
                         set_number, size_control, thread, penalty_size, name):
    """
    Function to launch the on the 3 exon list given by the user
    :param up_exon_list:  (list of exon instance) the list of exon up-regulated in a particular condition given
    by the user
    :param down_exon_list: (list of exon instance) the list of exon down-regulated in a particular condition given
    by the user
    :param up_and_down_list: (list of exon instance) the list of exon up and down regulated in a particular condition
    given by the user
    :param output_folder: (string) the path where your result will be created
    :param enrichment: (boolean) a boolean value indicating if the user want to do an enrichment analysis
    :param exon_type: (string) the type of the exon of the control sets
    :param set_number: (int)  the number of control sets to create
    :param size_control:  (boolean) a boolean value indicating if the user wants (True) or not (False) that the exons of
    the control sets share a similar length than his exons
    :param thread: (int) the number of processor the user wants to use to compute the control sets of exons
    :param penalty_size: (int) size in nucleotide of exon that will be considered as short when frequencies of codon,
    amino acid... are calculated and thus will have less weight than the others
    :param name: (string) corresponding to the name of the set of exon that will figure in the graphics produced by the
    program
    """
    os.mkdir(output_folder + "comparative_graphs/")
    codon_frequencies_up, aa_frequencies_up, nt_up, nt_exons_up = \
        main_program(up_exon_list, output_folder + "up/", enrichment, exon_type, set_number, size_control, thread,
                     penalty_size, "up " + name)
    codon_frequencies_down, aa_frequencies_down, nt_down, nt_exons_down = \
        main_program(down_exon_list, output_folder + "down/", enrichment, exon_type, set_number, size_control, thread,
                     penalty_size, "down " + name)

    if up_and_down_list is not None:
        codon_frequencies, aa_frequencies, nt_down, nt_exons = \
            main_program(up_and_down_list, output_folder + "up_and_down/", enrichment, exon_type, set_number,
                         size_control, thread, penalty_size, "up et down " + name)

    create_a_graphic_up_vs_down(codon_frequencies_up, codon_frequencies_down, "codons",
                                output_folder + "comparative_graphs/", name, nt_up, nt_down)
    create_a_graphic_up_vs_down(aa_frequencies_up, aa_frequencies_down, "acides amines",
                                output_folder + "comparative_graphs/", name)
    create_an_up_vs_down_comparison(nt_exons_up, nt_exons_down, output_folder + "comparative_graphs/")

    #up_down_graphic_maker_with_stats(up_exon_list, down_exon_list, output_folder + "comparative_graphs/", name)



def boxplot_nucleotid_rich_maker(output, up_list, down_list, up_down_list):
    """

    :param output: (string) the folder where the graphics will be created
    :param up_list: (list of exon instance) the list of exon up-regulated in a particular condition given
    by the user
    :param down_list: list of exon instance) the list of exon down-regulated in a particular condition given
    by the user
    :param up_down_list: (list of exon instance) the list of exon up and down regulated in a particular condition
    given by the user
    """
    os.mkdir(output + "nucleotides_information_figures/")
    output += "nucleotides_information_figures/"
    up_dict = up_list.codon_rich_exon_analyser()
    if down_list is not None:
        down_dict = down_list.codon_rich_exon_analyser()
    else:
        down_dict = None
    if up_down_list is not None:
        up_down_dict = up_down_list.codon_rich_exon_analyser()
    else:
        up_down_dict = None
    boxplot_nucleotide_rich_or_poor_codon(output, up_dict, down_dict, up_down_dict, "rich")
    boxplot_dinucleotide_rich_or_poor_codon(output, up_dict, down_dict, up_down_dict, "rich")

    up_dict = up_list.codon_poor_exon_analyser()
    if down_list is not None:
        down_dict = down_list.codon_poor_exon_analyser()
    else:
        down_dict = None
    if up_down_list is not None:
        up_down_dict = up_down_list.codon_poor_exon_analyser()
    else:
        up_down_dict = None

    boxplot_nucleotide_rich_or_poor_codon(output, up_dict, down_dict, up_down_dict, "poor")

    up_dict = up_list.codon_plus_exon_analyser()
    if down_list is not None:
        down_dict = down_list.codon_plus_exon_analyser()
    else:
        down_dict = None
    if up_down_list is not None:
        up_down_dict = up_down_list.codon_plus_exon_analyser()
    else:
        up_down_dict = None

    boxplot_dinucleotide_rich_or_poor_codon(output, up_dict, down_dict, up_down_dict, "plus")

    up_dict = up_list.codon_last_nt_proportion()
    if down_list is not None:
        down_dict = down_list.codon_last_nt_proportion()
    else:
        down_dict = None
    if up_down_list is not None:
        up_down_dict = up_down_list.codon_last_nt_proportion()
    else:
        up_down_dict = None

    boxplot_nucleotide_rich_or_poor_codon(output, up_dict, down_dict, up_down_dict, "last_nt")


def making_p_value_dictionary(dic_up_codon, dic_down_codon, dic_up_aa, dic_down_aa):
    """

    :param dic_up_codon: (dictionary) for each possible codon gives its frequency for all up-regulated exons
    :param dic_down_codon: (dictionary) for each possible codon gives its frequency for all down-regulated exons
    :param dic_up_aa: (dictionary) for each possible aa gives its frequency for all up-regulated exons
    :param dic_down_aa: (dictionary) for each possible aa gives its frequency for all down-regulated exons
    :return: 6 dictionaries of float:
    dic_p_value_codon : gives the p-value obtained by comparing the frequencies of a codon between the up and down
    regulated set of exon (for all possible codon)
    dic_p_adjust_codon : gives the corrected p-value (BH method ) obtained by comparing the frequencies of a codon
    between the up and down regulated set of exon (for all possible codon)
    dic_p_value_aa : gives the p-value obtained by comparing the frequencies of a aa between the up and down
    regulated set of exon (for all possible aa)
    dic_p_adjust_aa :gives the corrected p-value (BH method ) obtained by comparing the frequencies of a aa
    between the up and down regulated set of exon (for all possible aa)
    """
    rstats = importr('stats')
    # calculation of codon p-value
    p_value_codon = []
    parameter_difference_codon = []
    count = 1
    file_dir = os.path.dirname(os.path.realpath(__file__))
    for codon in dic_up_codon.keys():
        print_progress(count, 64, prefix='Statistical comparison up_vs_down codon', suffix='', decimals=1, bar_length=100)
        count += 1
        str_up = str(dic_up_codon[codon]).replace("[", "").replace("]", "").replace(" ", "")
        str_down = str(dic_down_codon[codon]).replace("[","").replace("]","").replace(" ","")
        res = commands.getoutput("Rscript --vanilla " + file_dir + "/comparison_up_and_down_freq.R " + str_up + " " + str_down)
        if "converged" in res:
            p_value_codon.append("NaN")
            parameter_difference_codon.append(["NaN", "NaN", "NaN", "NaN", "NaN", "NaN"])
        try:
            p_value_codon.append(float(res.split("p-value= ")[1].replace("\n","").strip().split(" ")[0]))
            tmp = res.split("p-value= ")[1].strip()
            sigma_up = tmp.split("=")[1].split("\n")[0]
            mu_up = tmp.split("=")[2].split("\n")[0]
            nu_up = tmp.split("=")[3].split("\n")[0]
            sigma_down = tmp.split("=")[4].split("\n")[0]
            mu_down = tmp.split("=")[5].split("\n")[0]
            nu_down = tmp.split("=")[6].split("\n")[0]
            parameter_difference_codon.append([sigma_up, mu_up, nu_up, sigma_down, mu_down, nu_down])
        except:
            p_value_codon.append("NaN")
            parameter_difference_codon.append(["NaN", "NaN", "NaN", "NaN", "NaN", "NaN"])
    p_adjust_codon = rstats.p_adjust(FloatVector(p_value_codon), method="BH")

    # calculation of aa_p_value
    p_value_aa = []
    parameter_difference_aa = []
    count = 1
    for aa in dic_up_aa.keys():
        print_progress(count, 21, prefix='Statistical comparison up_vs_down aa', suffix='', decimals=1, bar_length=100)
        count += 1
        str_up = str(dic_up_aa[aa]).replace("[", "").replace("]", "").replace(" ", "")
        str_down = str(dic_down_aa[aa]).replace("[","").replace("]","").replace(" ","")
        res = commands.getoutput("Rscript --vanilla " + file_dir + "/comparison_up_and_down_freq.R " + str_up + " " + str_down)
        try:
            p_value_aa.append(float(res.split("p-value= ")[1].strip().split(" ")[0]))
            tmp = res.split("p-value= ")[1].strip()
            sigma_up = tmp.split("=")[1].split("\n")[0]
            mu_up = tmp.split("=")[2].split("\n")[0]
            nu_up = tmp.split("=")[3].split("\n")[0]
            sigma_down = tmp.split("=")[4].split("\n")[0]
            mu_down = tmp.split("=")[5].split("\n")[0]
            nu_down = tmp.split("=")[6].split("\n")[0]
            parameter_difference_aa.append([sigma_up, mu_up, nu_up, sigma_down, mu_down, nu_down])
        except:
            p_value_aa.append("NaN")
            parameter_difference_aa.append(["NaN", "NaN", "NaN", "NaN", "NaN", "NaN"])
    p_adjust_aa = rstats.p_adjust(FloatVector(p_value_aa), method="BH")

    # turn the list into dictionaries for codon
    dic_p_value_codon = {}
    dic_p_adjust_codon = {}
    dic_param_codon = {}
    count = 0
    print len(parameter_difference_codon)
    for codon in dic_up_codon.keys():
        dic_p_value_codon[codon] = p_value_codon[count]
        dic_p_adjust_codon[codon] = p_adjust_codon[count]
        dic_param_codon[codon] = parameter_difference_codon[count]
        count += 1

    # turn the list into dictionaries for amino acid
    dic_p_value_aa = {}
    dic_p_adjust_aa = {}
    dic_param_aa = {}
    count = 0
    for aa in dic_up_aa.keys():
        dic_p_value_aa[aa] = p_value_aa[count]
        dic_p_adjust_aa[aa] = p_adjust_aa[count]
        dic_param_aa[aa] = parameter_difference_aa[count]
        count += 1

    return dic_p_value_codon, dic_p_adjust_codon, dic_p_value_aa, dic_p_adjust_aa, dic_param_codon, dic_param_aa

def up_down_graphic_maker_with_stats(up_exon_list, down_exon_list, output, name=None):
    if name is None:
        name = "du set d'interet"
    dic_up_codon = up_exon_list.codon_frequency_per_exon()
    dic_down_codon = down_exon_list.codon_frequency_per_exon()
    dic_up_aa = up_exon_list.amino_acid_frequency_per_exon()
    dic_down_aa = down_exon_list.amino_acid_frequency_per_exon()
    dic_p_value_codon, dic_p_adjust_codon, dic_p_value_aa, dic_p_adjust_aa, dic_param_codon, dic_param_aa = \
        making_p_value_dictionary(dic_up_codon, dic_down_codon, dic_up_aa, dic_down_aa)

    graphic_up_vs_down_stat(dic_up_codon, dic_down_codon, "codons", dic_p_value_codon, dic_p_adjust_codon, output, name)
    graphic_up_vs_down_stat(dic_up_aa, dic_down_aa, "acides amines", dic_p_value_aa, dic_p_adjust_aa, output, name)
    graphic_up_vs_down_param(dic_param_codon, "codons", output)
    graphic_up_vs_down_param(dic_param_aa, "aa", output)

    fichier = open(output + "res_stat.txt", "w")
    fichier.write("dic_up_codon = " + str(dic_up_codon) + "\n\n")
    fichier.write("dic_down_codon = " + str(dic_down_codon) + "\n\n")
    fichier.write("dic_p_value_codon = " + str(dic_p_value_codon) + "\n\n")
    fichier.write("dic_p_adjust_codon = " + str(dic_p_adjust_codon) + "\n\n")
    fichier.write("dic_up_aa = " + str(dic_up_aa) + "\n\n")
    fichier.write("dic_down_aa = " + str(dic_down_aa) + "\n\n")
    fichier.write("dic_p_value_aa = " + str(dic_p_value_aa) + "\n\n")
    fichier.write("dic_p_adjust_aa = " + str(dic_p_adjust_aa) + "\n\n")
    fichier.close()

    fichier2 = open(output + "res_parameters.txt", "w")
    fichier2.write("feature;sigma_up;mu_up;nu_up;sigma_down;mu_down;nu_down\n")
    for codon in dic_param_codon:
        fichier2.write(codon + ";" +  str(dic_param_codon[codon]).replace("[","").replace("]","").replace(",",";").replace(".",",").replace("'","") + "\n")
    for aa in dic_param_aa:
        fichier2.write(aa + ";" + str(dic_param_aa[aa]).replace("[","").replace("]","").replace(",",";").replace(".",",").replace("'","") + "\n")
    fichier2.close()



def launcher():
    """
    function that contains a parser to launch the program
    """
    # description on how to use the program
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description="""From a request file with the following shape :
    exon_name_1  chr_number  chromosomal_coordinate_1 chromosomal_coordinate_2
    exon_name_2  chr_number  chromosomal_coordinate_1 chromosomal_coordinate_2

With many exons as you want, gave the genomic sequence, the cds and the
codon/anticodon/amino_acid composition of those exons. Values of this files
have to be tab separated.

It is possible to specify or not an output folder. If an output folder is specify,
all the result files will be created in it. Otherwise, an output folder named
Exon_analysis is created in your current directory.

If you want to test whether your exons set is enriched (or impoverished) in particular codons, amino acids or in amino
acid nature, just precise it with the argument --enrichment True.
Moreover you have the possibility to chose the size of the control sets that will be use to test the enrichment or
impoverishment of your interest set. The bigger the control sets is, the more accurate the results will be.
    """)
    # Arguments for the parser

    parser.add_argument('--input', dest='input',
                        help="""your request containing an exon name, its chromosome number and its
                        chromosomal coordinates""", default=None)
    parser.add_argument('--up', dest='up', help="""your request file containing a list of exons up-regulated in a
                        certain condition it must be followed by the argument --down and a file containing the list of
                        the exons down-regulated in the same condition""")
    parser.add_argument('--down', dest='down', help="""your request file containing a list of exons down-regulated in a
                        certain condition it must be followed by the argument --up and a file containing the list of the
                        exons up-regulated in the same condition""")
    parser.add_argument('--up_and_down', dest='up_and_down', help="""if you have already set the up and down parameters
    you may want to analysis the up and down set together : in that case, type : --up_and_down yes """, default="no")

    parser.add_argument('--output', dest='output', help="An output folder",
                        default=".")

    parser.add_argument('--enrichment', dest='enrichment', help='perform the enrichment analysis', default=False)

    parser.add_argument('--exon_type', dest='exon_type', help='exon type of the control set', default=False)

    parser.add_argument('--set_number', dest='set_number', help='the number of control sets to be created',
                        default=10000)

    parser.add_argument('--size_control', dest='size_control',
                        help='enable control size : the exon of the control sets will have approximately the same size '
                             'than those given by the user',
                        default=True)

    parser.add_argument('--thread', dest='thread',
                        help='number of thread you want to use to create the control data sets',
                        default=1)

    parser.add_argument('--complete_exon', dest='complete_exon',
                        help='True if you want to analyse only full exon, and not potential overlapping exon',
                        default=True)

    parser.add_argument('--penalty_size', dest='penalty',
                        help="""size in nucleotides of exons that will be considered as short when frequencies of
                        codons, amino acid... are calculated and thus will have less weight than the others""",
                        default=51)

    parser.add_argument("--correspondence", dest="correspondence",
                        help="""If the exon name (first name of the input column) corresponds to a gene name following
                               by the exon position of this gene, this argument allow to check if the exon name
                               corresponds to the gene name and the exon position (useful in case of different
                               annotations)""",
                        default=False)

    parser.add_argument("--duplication", dest="duplication",
                        help="""True if we want to allow duplication interval, False else""",
                        default=False)

    parser.add_argument("--set_name", dest="set_name",
                        help="""the name of the set that you want in you figures""",
                        default="du set d'interet")

    parser.add_argument("--alt_figures", dest='alt_figures',
                        help="yes if you want to create alternative figures for your analysis", default="no")
    parser.add_argument("--summary", dest='summary', help="yes if you want to create a summary for your figure",
                        default="no")

    args = parser.parse_args()  # parsing arguments
    try:
        args.penalty = int(args.penalty)
    except TypeError:
        print "Penalty_size should be a number"
        print "Terminating..."
        exit()

    if args.exon_type is False and args.enrichment is True:
        print "To process enrichment analysis an exon_type have to be given to the program"
        print "Terminating..."
        exit()

    try:
        args.thread = int(args.thread)
    except TypeError:
        print "Invalid number of thread"
        print "Terminating..."
        exit()

    try:
        args.set_number = int(args.set_number)
    except TypeError:
        print "Invalid set_number"
        print "Terminating..."
        exit()

    if args.input is not None:
        if os.path.isfile(args.input) is False:
            print "The input file doesn't exist"
            print "Terminating..."
            exit()
    else:
        if not os.path.isfile(args.up) or not os.path.isfile(args.down):
            print "The up or/and the down exon file doesnt exist"
            print "Terminating..."
            exit()

    if args.up_and_down != "yes" and args.up_and_down != "no":
        print "Wrong value for \"up_and_down\" parameter"
        print "Setting the parameter up_and_down to : no"
        args.up_and_down = "no"

    if args.up_and_down == "yes" and args.input is not None:
        print "You can't concatenate the up-regulated and the down-regulated exons of your project because you only " \
              "gave one file to the program"
        print "Setting the parameter up_and_down to : no"
        args.up_and_down = "no"

    args.output = checking_output_path(args.output)

    if args.complete_exon == "True":
        args.complete_exon = True

    if args.correspondence == "True":
        args.correspondence = True

    if args.enrichment == "False":
        args.enrichment = False

    if args.enrichment == "True":
        args.enrichment = True

    if args.duplication == "False":
        args.duplication = False

    if args.duplication == "True":
        args.duplication = True

    if args.alt_figures != "yes" and args.alt_figures != "no":
        print "Wrong value for the parameter 'alt_figure' "
        print "Setting its value to 'no'"

    if args.summary != "yes" and args.summary != "no":
        print "Wrong value for the parameter 'summary' "
        print "Setting its value to 'no'"

    if len(commands.getoutput("ls " + args.output).split("\n")) > 1:
        rep = raw_input("the folder " + str(args.output) + " already contains result, do you want to remove them ?"
                                                           "(y/n) : ")
        while rep != "y" and rep != "n":
            rep = raw_input("the folder " + str(args.output) + " already contains result, do you want to remove them?"
                                                               "(y/n) : ")
        if rep == "y":
            commands.getoutput("rm -rf " + str(args.output) + "*")
        else:
            print "exiting..."
            exit()
    if args.input is not None:
        full_exon_list, input_content = execute_request(args.input, args.output, args.complete_exon,
                                                        args.correspondence, args.duplication)

        if args.enrichment is True:
            main_program(full_exon_list, args.output, args.enrichment, args.exon_type, args.set_number,
                         args.size_control, args.thread, args.penalty, args.set_name)
            boxplot_nucleotid_rich_maker(args.output, full_exon_list, None, None)
        if args.summary == "yes":
            creating_composition_summary(args.output, full_exon_list, None, None, args.penalty)

        if args.alt_figures == "yes":
            creating_recap_graphics(args.output, full_exon_list, None, None)
    else:
        if args.up_and_down == "no":
            os.mkdir(args.output+"up/")
            os.mkdir(args.output + "down/")
            # executing the program
            up_exon_list, input_content = execute_request(args.up, args.output + "up/", args.complete_exon,
                                                          args.correspondence, args.duplication)
            down_exon_list, input_content = execute_request(args.down, args.output + "down/", args.complete_exon,
                                                            args.correspondence, args.duplication)


            if args.enrichment is True:
                # executing the program
                up_and_down_launcher(up_exon_list, down_exon_list, None, args.output, args.enrichment, args.exon_type,
                                     args.set_number, args.size_control, args.thread, args.penalty, args.set_name)

            if args.summary == "yes":
                creating_composition_summary(args.output, up_exon_list, down_exon_list, None,
                                             args.penalty)

            boxplot_nucleotid_rich_maker(args.output, up_exon_list, down_exon_list, None)
            if args.alt_figures == "yes":
                creating_recap_graphics(args.output, up_exon_list, down_exon_list, None)

        else:
            os.mkdir(args.output+"up/")
            os.mkdir(args.output + "down/")
            os.mkdir(args.output + "up_and_down/")
            # executing the program
            up_exon_list, input_content_up = execute_request(args.up, args.output+"up/", args.complete_exon,
                                                             args.correspondence, args.duplication)
            down_exon_list, input_content_down = execute_request(args.down, args.output + "down/", args.complete_exon,
                                                                 args.correspondence, args.duplication)
            all_content = input_content_up + input_content_down
            all_exon_list = concatenate_up_and_down_exon_list(args.output + "up_and_down/", up_exon_list,
                                                              down_exon_list, all_content)


            if args.enrichment is True:
                # executing the program
                up_and_down_launcher(up_exon_list, down_exon_list, all_exon_list, args.output, args.enrichment,
                                     args.exon_type, args.set_number, args.size_control, args.thread, args.penalty,
                                     args.set_name)
                boxplot_nucleotid_rich_maker(args.output, up_exon_list, down_exon_list, all_exon_list)

            if args.summary == "yes":
                creating_composition_summary(args.output, up_exon_list, down_exon_list, all_exon_list,
                                             args.penalty)

            if args.alt_figures == "yes":
                creating_recap_graphics(args.output, up_exon_list, down_exon_list, all_exon_list)


if __name__ == "__main__":
    launcher()  # launches the program
