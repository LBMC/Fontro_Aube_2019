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
from list_of_exon import *
sys.path.insert(0, os.path.realpath(os.path.dirname(__file__)).replace("/frequences", ""))
import config
from exon_class import *


##########################
# Functions
##########################

def connection():
    """
    :return: an object that contains all the information you need to connect to fasterDB
    """
    cnx = mysql.connector.connect(user=config.user, password=config.password, host=config.host,
                                  database=config.database, buffered=True)
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


def stop_message():
    """
    If the sanitization check (see the function below) detects a problem in the input fill, this function is called
    to inform the user that the program will stop
    """
    print "Input did not pass sanitization check."
    print "So die..."
    exit()


def sanitization_check(input_file):
    """
    :param input_file: the input fil containing the exons that the user wants to study
    If an input line is not ok, nothing is done, the program stops
    :return: the number of line of the request
    """
    request = open(input_file, "r")
    line = request.readline()
    counter = 1
    while line:
        line = line.split("\t")

        if len(line) != 4:
            stop_message()
        if str(line[1]) not in ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13",
                                "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y"]:
            stop_message()
        try:
            int(line[2])
            int(line[3])
        except ValueError:
            print "Input did not pass sanitization check."
            print "So die..."
            exit()
        line = request.readline()
        counter += 1
    request.close()
    return counter


def copy_input_name(input_file, path):
    """
    :param input_file: the name of the input file
    :param path: the path were the copy of the input will be created
    :return: a tuple containing the name of the copy of the input file that will be drop in the output directory
    The new name is calculated to make sur this new file will not erase an existing one
    """
    file_name = input_file + "_copy_1.csv"
    counter = 1
    while os.path.isfile(path + file_name):
        counter += 1
        file_name = input_file + "_copy_" + str(counter) + ".txt"
    return counter, file_name


def print_progress(iteration, total, prefix='', suffix='', decimals=1, barlength=100):
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
    filled_length = int(round(barlength * iteration / float(total)))
    bar = "#" * filled_length + '-' * (barlength - filled_length)
    sys.stdout.write('\r%s |%s| %s%s %s' % (prefix, bar, percent, '%', suffix)),
    if iteration == total:
        sys.stdout.write('\n')
    sys.stdout.flush()


def execute_request(input_file, name, same_gene, complete_exon, count, correspondence, output):
    """
    Calculate for each interval (that normally corresponds to an exon) all exons \
    that are in this interval, their sequences (genomic, cds peptide) and their composition (in codons, \
    corresponding anticodons and amino acids)

    :param input_file: (string) the file containing the input of the exon
    :param name: (string) the name of the file thta will be created
    :param same_gene: (boolean) True if the input file contains exon comming from the same gene, false else
    :param complete_exon: (boolean) True if you want to retrieve only exon matching exactly with the intervals \
    given in ``input_file``
    :param count: (boolean) True if you want to write count amino acid and count files, false if you want \
    to write frequency file of amino acids and codons.
    :param correspondence: a boolean value indicating if the user wants the gene name and the exon position on this gene
    :param output: (string) folder were the file will be created
    """
    input_size = sanitization_check(input_file)  # checks if the input is correct
    request = open(input_file, "r")
    line = request.readline()
    cnx = connection()  # loading all the information to connect to fasterDB
    full_exon_list = ListExon()
    counter = 0
    while line:  # reading the input file
        if input_size > 300:
            print_progress(counter, input_size)
            counter += 1
        line = line.split("\t")
        line[2] = int(line[2])
        line[3] = int(line[3])
        list_of_exon = retrieve_matching_exon(cnx, line[1], line[2:4], line[0])
        if list_of_exon is not None:  # if the interval given by the user match exons we :
            for i in range(len(list_of_exon.exon_list)):
                list_of_exon.exon_list[i].retrieve_gene_name(cnx)  # retrieve the name of the genes that contains it
                list_of_exon.exon_list[i].retrieve_exon_sequences(cnx)  # retrieve its sequences
                list_of_exon.exon_list[i].calculate_exon_coverage_on_input_and_input_coverage_on_exon(line[2:4])
                # resize the sequences if the user only select a part of the exon
                list_of_exon.exon_list[i].resize_cds_sequence(line[2:4])
                # compute all the exon, corresponding anticodons, amino acids within the exon
                list_of_exon.exon_list[i].found_codon_anticodon_amino_acid_and_aa_nature()
            # we increment those new exon with no lacking information in a list containing all exons
            full_exon_list.exon_list += list_of_exon.exon_list
        line = request.readline()
    request.close()
    # writing output files
    if complete_exon is True and correspondence is False:
        i = 0
        while i < len(full_exon_list.exon_list):
            if full_exon_list.exon_list[i].matching[1] != 100.0 or full_exon_list.exon_list[i].matching[2] != 100.0:
                del(full_exon_list.exon_list[i])
                i -= 1
            i += 1
    elif complete_exon is True and correspondence is True:
        i = 0
        while i < len(full_exon_list.exon_list):
            if full_exon_list.exon_list[i].matching[1] != 100.0 or full_exon_list.exon_list[i].matching[2] != 100.0 \
                or full_exon_list.exon_list[i].exon_name.split("_")[0] != full_exon_list.exon_list[i].gene_name\
                    or full_exon_list.exon_list[i].exon_name.split("_")[1] != full_exon_list.exon_list[i].exon_number:
                    del (full_exon_list.exon_list[i])
                    i -= 1
            i += 1
    print "list len = " + str(len(full_exon_list.exon_list))
    if same_gene == "True":
        if full_exon_list.exon_gene() is False:
            print "Warning : your exon set do not come from the same gene !"
        full_exon_list.write_codon_frequency_calculator(name, count, output)
        full_exon_list.write_aa_frequency_calculator(name, count, output)
    else:
        if full_exon_list.exon_gene() is True:
            print "Warning : your exon set come from a only gene, " \
                  "the parameter same_gene True is more appropriate here !"
        full_exon_list.write_codon_weighted_frequency_calculator(17, name, output)
        full_exon_list.write_aa_weighted_frequency_calculator(17, name, output)


def main():
    """
    function that contains a parser to launch the program
    """
    # description on how to use the program
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description="""From a request file with the following shape :
    exon_name_1  chr_number  chromosomal_coordinate_1 chromosomal_coordinate_2
    exon_name_2  chr_number  chromosomal_coordinate_1 chromosomal_coordinate_2

With many exons as you want.
######## It give for each codon their frequency in the exon set ###############
######## If the exon set do not come from the same gene you should precise it ##########

It is possible to specify or not an output folder. If an output folder is specify,
 all the result files will be created in it. Otherwise, an output folder named
 Exon_analysis is created in your current directory
    """,
                                     usage='%(prog)s --input input_file.txt [--output an output folder] ')
    # Arguments for the parser
    required_name = parser.add_argument_group('required arguments')
    required_name.add_argument('--input', dest='input', required=True,
                               help="""your request containing an exon name, its chromosome number and its
                                chromosomal coordinates""")

    parser.add_argument('--name_file', dest='name_file', help="A name for your file",
                        default="result")

    parser.add_argument('--same_gene', dest='same_gene',
                        help="True/False according to the fact that the input exons come from "
                             "the same or a different gene",
                        default="True")

    parser.add_argument('--complete_exon', dest='complete_exon',
                        help="True is you want to analyse only complete exon",
                        default=True)

    parser.add_argument('--count', dest='count',
                        help="True if you want to have the count for each codon instead of the frequency",
                        default=False)

    parser.add_argument('--output', dest='output', help="path were the frequency files will be created",
                        default="result")

    args = parser.parse_args()  # parsing arguments
    execute_request(args.input, args.name_file, args.same_gene, args.complete_exon, args.count, True, args.output)


if __name__ == "__main__":
    main()  # launches the program
