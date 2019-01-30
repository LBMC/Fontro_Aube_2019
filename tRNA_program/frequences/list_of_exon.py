import sys
import os
sys.path.insert(0, os.path.realpath(os.path.dirname(__file__)).replace("/frequences", ""))
from dictionnary import codon2aminoAcid  # Links each codon to its respective anticodon
from dictionnary import amino_acid2codon


class ListExon:
    """
    A class that contains a list of exon.
    Designed to wrap a set of exon and to links it with method
    """

    def __init__(self):
        """
        Initialisation of the attribute exon_list that contains a list of exon
        """
        self.exon_list = list()

    def write_codon_frequency_calculator(self, name, count, output="."):
        """
        Calculate the frequency of the codons in the gene and write those information in a file
        The exon have to be a part of the same gene
        """
        file_name = output + "/codon_frequency_" + str(name) + ".txt"
        codon_frequency_file = open(file_name, "w")
        codon_frequency_file.write("Codon\tFrequency\n")
        codon_list = list()
        for i in range(len(self.exon_list)):
            if len(self.exon_list[i].codon) > 0:
                codon_list += self.exon_list[i].codon
        for key in codon2aminoAcid.keys():
            if count is False:
                codon_frequency_file.write(str(key)+"\t"+str(float(codon_list.count(key)) / len(codon_list))+"\n")
            else:
                codon_frequency_file.write(str(key) + "\t" + str(float(codon_list.count(key))) + "\n")
        codon_frequency_file.close()

    def write_aa_frequency_calculator(self, name, count, output="."):
        """
        Calculate the frequency of the amino acid in the gene and write those information in a file
        The exon have to be a part of the same gene
        """
        file_name = output + "/aa_frequency_" + str(name) + ".txt"
        aa_frequency_file = open(file_name, "w")
        aa_frequency_file.write("AA\tFrequency\n")
        aa_list = list()
        # fusion of the exon
        for i in range(len(self.exon_list)):
            if len(self.exon_list[i].amino_acid) > 0:
                aa_list += self.exon_list[i].amino_acid
        for key in amino_acid2codon.keys():
            if count is False:
                aa_frequency_file.write(str(key) + "\t" + str(float(aa_list.count(key)) / len(aa_list)) + "\n")
            else:
                aa_frequency_file.write(str(key) + "\t" + str(float(aa_list.count(key))) + "\n")
        aa_frequency_file.close()

    def write_codon_weighted_frequency_calculator(self, length_penalty, name, output="."):
        """
        weight the frequency of a codon by the length of the exon and by the number of codon
        for instance the set of codon : ATGCCA  ATGATG  GAGCCA CTG
        will have an ATG frequency of ATG = (0.5 + 1 + 0 + 0)/4 = 0.375
        """
        file_name = output + "/codon_frequency_" + str(name) + ".txt"
        codon_frequency_file = open(file_name, "w")
        codon_frequency_file.write("Codon\tFrequency\n")
        dic = {}
        for key in codon2aminoAcid.keys():
            dic[key] = 0.
            c = 0.
            for i in range(len(self.exon_list)):
                if len(self.exon_list[i].codon) > 0:
                    if len(self.exon_list[i].codon) > length_penalty - 1:
                        c += 1.
                        dic[key] += float(self.exon_list[i].codon.count(key)) / len(self.exon_list[i].codon)
                    else:
                        c += float(len(self.exon_list[i].codon)) / length_penalty
                        dic[key] += float(self.exon_list[i].codon.count(key)) / \
                            len(self.exon_list[i].codon) * len(self.exon_list[i].codon) / length_penalty
            dic[key] = dic[key]/c
            codon_frequency_file.write(str(key) + "\t" + str(dic[key]) + "\n")
        codon_frequency_file.close()

    def write_aa_weighted_frequency_calculator(self, length_penalty, name, output="."):
        """
        weight the frequency of a codon by the length of the exon and by the number of codon
        for instance the set of codon : ATGCCA  ATGATG  GAGCCA CTG
        will have an ATG frequency of ATG = (0.5 + 1 + 0 + 0)/4 = 0.375
        """
        file_name = output + "/aa_frequency_" + str(name) + ".txt"
        aa_frequency_file = open(file_name, "w")
        aa_frequency_file.write("Codon\tFrequency\n")
        dic = dict()
        for key in amino_acid2codon.keys():
            dic[key] = 0.
            c = 0.
            for i in range(len(self.exon_list)):
                if len(self.exon_list[i].amino_acid) > 0:
                    if len(self.exon_list[i].amino_acid) > length_penalty - 1:
                        c += 1.
                        dic[key] += float(self.exon_list[i].amino_acid.count(key)) / len(self.exon_list[i].amino_acid)
                    else:
                        c += float(len(self.exon_list[i].amino_acid)) / length_penalty
                        dic[key] += float(self.exon_list[i].amino_acid.count(key)) / \
                            len(self.exon_list[i].amino_acid) * len(self.exon_list[i].amino_acid) / length_penalty
            dic[key] = dic[key]/c
            aa_frequency_file.write(str(key) + "\t" + str(dic[key]) + "\n")
        aa_frequency_file.close()

    def exon_gene(self):
        """
        :return: true if the set of exon are coming from the same gene, else return false
        """
        gene = self.exon_list[0].gene_name
        for i in range(1, len(self.exon_list)):
            if gene != self.exon_list[i].gene_name:
                return False
        return True
