#!/usr/bin/env python2

# Import
import frequency_finder
import os


# functions
def retrieve_control_exon(cnx, exon_type):
    """
    :param cnx: (mysql.connector instance) the information necessary to connect to fasterDB
    :param exon_type: (string) the type of control exon analyzed
    :return: (string) the name of the file created
    """
    outfolder = os.path.realpath(os.path.dirname(__file__)) + "/input"
    outfile = "%s/input_%s.txt" % (outfolder, exon_type)
    if not os.path.isdir(outfolder):
        os.mkdir(outfolder)
    cursor = cnx.cursor()
    query = ("""
             SELECT t2.official_symbol, t1.pos_sur_gene, t1.chromosome, t1.start_sur_chromosome, t1.end_sur_chromosome 
             FROM hsapiens_exonsstatus_improved t1, IMPORT_FasterDB_genes t2
             WHERE t1.id_gene = t2.id
             AND t1.exon_types LIKE \"%{}%\"""".format(exon_type))
    cursor.execute(query)
    with open(outfile, "w") as my_file:
        for exon in cursor:
            my_file.write("%s_%s\t%s\t%s\t%s\n" % (exon[0], exon[1], exon[2], exon[3], exon[4]))
    return outfile


def run_frequency_finder(input_file, exon_type):
    """
    Run frequency finder.

    :param input_file: (string) a file containing exons
    :param exon_type: (string) the type of exon of interest
    :return:
    """
    outfolder = os.path.realpath(os.path.dirname(__file__)) + "/frequence"
    if not os.path.isdir(outfolder):
        os.mkdir(outfolder)
    frequency_finder.execute_request(input_file, exon_type, False, True, False, True, outfolder)


def main():
    exon_types = ["ACE", "CCE", "FCE", "LCE"]
    cnx = frequency_finder.connection()
    for exon_type in exon_types:
        input_file = retrieve_control_exon(cnx, exon_type)
        run_frequency_finder(input_file, exon_type)


if __name__ == "__main__":
    main()
