#!/usr/bin/python3


"""
Description:

The goal of this script is to:
    1. Create a lot of random fasta enriched in a particular unit (*nucleotide*, *di-nucleotide* or *feature*) and calculates the mean frequency of this unit in the fasta file
    2. Create a lot of random fasta impoverished in **the same unit** (*nucleotide*, *di-nucleotide* or *feature*) and calculates the mean frequency of this unit in the fasta file

And then compare their frequency for this unit for each couple of random fasta enriched and impoverished \
(for this unit). To compare their frequency, the relative frequency is computed

"""

import subprocess
from Bio import SeqIO
import numpy as np
import argparse
import os
import sys
import stretch_evalutator
import config


#stretches = [[4, 5], [5, 6], [6, 7], [7, 8], [8, 9], [9, 10], [10, 11], [11, 12]]
stretches = [[4, 5], [5, 6], [6, 7], [7, 8], [8, 9], [9, 10], [10, 11], [11, 12]]

feature_dic = {
    "Very-small": ["A", "C", "G", "S"],
    "Small#2": ["A", "C", "D", "G", "N", "P", "S", "T"],
    "Large" : ["F", "I", "K", "L", "M", "R", "W", "Y"],
    "Disorder-promoting#1": ["A", "E", "G", "K", "P", "Q", "R", "S"],
    "Order-promoting#1": ["C", "F", "I", "L", "N", "W", "V", "Y"],
    "Disorder-promoting#2": ["A", "E", "G", "K", "P", "Q", "S"],
    "Order-promoting#2": ["C", "F", "H", "I", "L", "M", "N", "W", "V", "Y"],
    "Polar-uncharged#1": ["C", "N", "Q", "S", "T", "Y"],
    "Polar-uncharged#2": ["N", "Q", "S", "T", "Y"],
    "Charged#1": ["R", "H", "K", "D", "E"],
    "Charged#2": ["R", "K", "D", "E"],
    "Hydrophilic#1": ["D", "E", "K", "N", "Q", "R"],
    "Hydrophobic#1": ["A", "C", "F", "I", "L", "M", "V"],
    "Neutral": ["G", "H", "P", "S", "T", "Y"],
    "Hydroxylic": ["S", "T", "Y"],
    "Negatively-charged": ["D", "E"],
    "Positively-charged#1": ["R", "H", "K"],
    "Positively-charged#2": ["R", "K"]
}


def fasta_generator_dependant_feature(feature, freq, output):
    """
    Launch the genesis of sequence with dependant enrichment and impoverishment of \
    the 2 feature given in the list feature ``feature``.

    :param feature: (list of string) list of feature
    :param freq:  (list of float) list of prop
    :param output: (string) path where the fasta file will be created
    :return:
    """

    cmd = ["python2"]
    cmd += [config.fg2fre]
    cmd += ["--output", output]
    cmd += ["--feature1", feature[0]]
    cmd += ["--feature2", feature[1]]
    cmd += ["--prop1", freq[0]]
    cmd += ["--prop2", freq[1]]
    subprocess.check_output(cmd, stderr=subprocess.STDOUT)


def translator(sequence):
    """
    Translate a nucleotide sequence. The ORF of the sequence starts at the first nucleotide and codons can't \
    be truncated.

    :param sequence: (string a nucleotide sequence
    :return: (string an amino acid sequence
    """
    codon2aminoacid = dict(TTT="F", TTC="F", TTA="L", TTG="L", CTT="L", CTC="L", CTA="L", CTG="L", ATT="I", ATC="I",
                           ATA="I", ATG="M", GTT="V", GTC="V", GTA="V", GTG="V", TCT="S", TCC="S", TCA="S", TCG="S",
                           CCT="P", CCC="P", CCA="P", CCG="P", ACT="T", ACC="T", ACA="T", ACG="T", GCT="A", GCC="A",
                           GCA="A", GCG="A", TAT="Y", TAC="Y", TAA="*", TAG="*", CAT="H", CAC="H", CAA="Q", CAG="Q",
                           AAT="N", AAC="N", AAA="K", AAG="K", GAT="D", GAC="D", GAA="E", GAG="E", TGT="C", TGC="C",
                           TGA="*", TGG="W", CGT="R", CGC="R", CGA="R", CGG="R", AGT="S", AGC="S", AGA="R", AGG="R",
                           GGT="G", GGC="G", GGA="G", GGG="G")
    new_seq = []
    for i in range(0, len(sequence), 3):
        if i+3 <= len(sequence):
            new_seq.append(sequence[i:i+3])
    aa_seq = ""
    for codon in new_seq:
        aa_seq += codon2aminoacid[codon]
    return aa_seq


def get_feature_frequency(sequence, feature):
    """
    Get the feature frequency of a nucleotide sequence. It's ORF starts at the first nucleotide and codons can't \
    be truncated.

    :param sequence: (string) a nucleotide sequence
    :param feature: (string) the feature for which we want to compute the frequency in the sequence
    :return: (float) the frequency of the feature ``feature``
    """
    aa_seq = translator(sequence)
    count = 0.
    for aa in feature_dic[feature]:
            count += aa_seq.count(aa)
    return (float(count) / len(aa_seq)) * 100


def get_dinucleotide_frequency(sequence, dnt):
    """
    Get the di-nucleotide frequency of a sequence

    :param sequence: (string) a nucleotide sequence
    :param dnt: (string) a di-nucleotide for which we want to calculate the frequency
    :return: (float) the frequency of the di-nucleotide ``dnt`` in ``sequence``
    """
    count = 0.
    for i in range(len(sequence) - 1):
        if sequence[i:i+2] == dnt:
            count += 1
    return (count / (len(sequence) - 1)) * 100


def get_nucleotide_frequency(sequence, nt):
    """
    Get the nucleotide frequency of a nucleotide sequence

    :param sequence: (string) a nucleotide sequence
    :param nt: (string) a nucleotide for which we want to compute teh frequency
    :return: (float) the frequency of the nucleotide ``nt`` in ``sequence``
    """
    iupac = {'Y': ['C', 'T'], 'R': ['A', 'G'], 'W': ['A', 'T'], 'S': ['G', 'C'], 'K': ['T', 'G'], 'M': ['C', 'A'],
             'D': ['A', 'G', 'T'], 'V': ['A', 'C', 'G'], 'H': ['A', 'C', 'T'], 'B': ['C', 'G', 'T']}
    count = 0.
    if nt in ["A", "C", "G", "T"]:
        for i in range(len(sequence)):
            if sequence[i] == nt:
                count += 1
    else:
        for i in range(len(sequence)):
            if sequence[i] in iupac[nt]:
                count += 1
    return (count / len(sequence)) * 100


def get_mean_frequency_in_fasta_file(fasta_file, type_unit, unit):
    """
    Give the mean frequency of the ``type_unit`` ``unit`` in the all set of sequence \
    contains in the fasta files.

    :param fasta_file: (string) path to a fasta files
    :param type_unit: (string) the unit type we want to study : it can be *dnt* or *feature or *nt*
    :param unit: (string) the unit for which we want to calculate the median frequency
    :return: the mean frequency of ``unit`` of all the sequence in ``fasta_file``
    """
    mean_list = []
    while len(mean_list) == 0:
        mean_list = []
        for record in SeqIO.parse(fasta_file, "fasta"):
            seq = str(record.seq)
            if type_unit == "feature":
                mean_list.append(get_feature_frequency(seq, unit))
            if type_unit == "dnt":
                mean_list.append(get_dinucleotide_frequency(seq, unit))
            if type_unit == "nt":
                mean_list.append(get_nucleotide_frequency(seq, unit))
    return np.mean(mean_list)


def get_mean_frequency_of_multiple_fasta(type_unit, unit, freq, iteration, output, iunit_type, iunit):
    """
    Create ``iteration`` fasta files having a frequency of ``freq`` of ``unit`` and gets the frequency of \
    ``unit`` in all the fasta files generated

    :param type_unit: (string) the unit type we want to study : it can be *dnt* or *feature or *nt*
    :param unit: (string) the unit for which we want to calculate the median frequency
    :param freq: (float) the frequency of ``unit`` in the fasta file that will be produced
    :param iteration: (int) the number of fasta files we want to create
    :param output: (string) path where the fasta_file will be created
    :param itype_unit:  (string) the unit type we want to study, it can be *dnt* or *feature or *nt*
    :param iunit: the unit of ``itype_unit`` for which we want to calculate the mean frequency in fasta files
    :return: (list of ``itaration`` float value) list of the mean of ``unit`` in the ``iteration`` fasta \
    files generated.
    """
    list_mean_fasta = []
    list_mean_fasta2 = []
    sizes = []
    list_count_stretch = [[] for j in range(len(stretches))]
    ilist = [[]]
    iunit2 = iunit.split(",")
    iunit_type2 = iunit_type.split(",")
    if len(iunit2) > 1:
        ilist = [[] for i in range(len(iunit2))]

    filename_fasta = "CCE_" + unit + "_" + str(freq)
    for i in range(iteration):
        sys.stdout.write("\rprogress: " + str(i + 1) + " / " + str(iteration))
        if i == iteration - 1:
            sys.stdout.write("\n")
        sys.stdout.flush()
        unit2 = unit.split(",")
        freq2 = freq.split(",")
        fasta_generator_dependant_feature(unit2, freq2, output)
        list_mean_fasta.append(get_mean_frequency_in_fasta_file(output + filename_fasta + ".fasta", type_unit, unit2[0]))
        list_mean_fasta2.append(get_mean_frequency_in_fasta_file(output + filename_fasta + ".fasta", type_unit, unit2[1]))
        if iunit_type is not None:
            for i in range(len(iunit2)):
                ilist[i].append(get_mean_frequency_in_fasta_file(output + filename_fasta + ".fasta", iunit_type2[i], iunit2[i]))
            if iunit_type2[0] == "nt":
                value_stretch, size = stretch_evalutator.stretch_calculator(output + filename_fasta + ".fasta", True, iunit_type2[0], iunit2[0], output, stretches)
                sizes.append(size[0])
                for k in range(len(list_count_stretch)): list_count_stretch[k].append(value_stretch[k])
    subprocess.check_call(["rm", output + filename_fasta + ".fasta"], stderr=subprocess.STDOUT)
    if iunit_type is None:
        return [list_mean_fasta, list_mean_fasta2]
    else:
        if iunit_type2[0] == "nt":
            if sizes[0] != sum(sizes) / len(sizes):
                print("The " + iteration + " fasta created did not have the same number of sequence")
            res =  [list_mean_fasta, list_mean_fasta2] + ilist + [list_count_stretch, sizes]
            return res
        else:
            res = [list_mean_fasta, list_mean_fasta2] + ilist + [None, None]
            return res


def get_relative_freq_values(list_high, list_low):
    """
    Get the relative frequencies of fasta_files having up and low frequency for a particular unit.

    :param list_high: (list of float) list of mean frequencies in a particular unit in a list of fasta file enriched \
    for this unit
    :param list_low: (list of float) list of mean frequencies in a particular unit in a list of fasta file impoverished \
    for this unit
    :return: (list of float) relative frequencies for each value of list_high and list low.

    """
    relative_freq = []
    for i in range(len(list_high)):
        relative_freq.append(((list_high[i] - list_low[i]) / list_low[i]) * 100)
    return relative_freq


def write_full_tsv_file(unit_type, unit, freq_high, freq_low, list_high_freq, list_low_freq,
                   iteration, output, iunit_type, iunit):
    """
    Create a tsv file.

    :param unit_type: (string) the unit type we want to study : it can be *dnt* or *feature or *nt*
    :param unit: (string) the unit for which we want to calculate the median frequency
    :param freq_high: (float) the frequency of ``unit`` for the fasta file having an high content of the unit of \
    interest
    :param freq_low: (float) the frequency of ``unit`` for the fasta file having a low content of the unit of \
    interest
    :param list_high_freq: (list of float) list of mean frequencies in a particular unit in a list of fasta file enriched \
    for this unit
    :param list_low_freq: (list of float) list of mean frequencies in a particular unit in a list of fasta file impoverished \
    for this unit
    :param iteration: (int) the number of fasta files we want to create
    :param output: (string) path where the fasta_file will be created
    :param iunit_type:  (string) the unit type we want to study, it can be *dnt* or *feature or *nt*
    :param iunit: the unit of ``itype_unit`` for which we want to calculate the mean frequency in fasta files
    """
    unit2 = unit.split(",")
    iunit2 = iunit.split(",")
    fname = "mutated"
    freq_high2 = freq_high.split(",")
    freq_low2 = freq_low.split(",")
    relative_freq = [get_relative_freq_values(list_high_freq[i], list_low_freq[i]) for i in range(len(list_high_freq))]
    mean_rel_freq = [np.mean(relative_freq[i]) for i in range(len(relative_freq))]
    std_rel_freq = [np.std(relative_freq[i], ddof=1)for i in range(len(relative_freq))]

    mean_high = [np.mean(list_high_freq[i]) for i in range(len(list_high_freq))]
    std_high = [np.std(list_high_freq[i], ddof=1) for i in range(len(list_high_freq))]

    mean_low = [np.mean(list_low_freq[i]) for i in range(len(list_high_freq))]
    std_low = [np.std(list_low_freq[i], ddof=1) for i in range(len(list_high_freq))]
    unit = unit.replace("_", "-")


    # T test between the frequency of the same feature in high and low fastas
    pval = [stretch_evalutator.r_ttest(list_high_freq[i], list_low_freq[i]) for i in range(len(list_high_freq))]

    with open(output + iunit_type + "_" + iunit + "_frequency_comparison_between_" + str(iteration) + "_" + fname + "_fasta_file-" + str(unit) + "_high:" + str(freq_high) + "_low:" + str(freq_low) + ".tsv", "w") as outfile:
        outfile.write("frequency_fasta" + fname + "_" + str(unit2[0]) + ":" + str(freq_high2[0]) + "\t")
        outfile.write("frequency_fasta" + fname + "_" + str(unit2[0]) + ":" + str(freq_low2[0]) + "\t")
        outfile.write("relative_frequency\t")
        outfile.write("frequency_fasta" + fname + "_" + str(unit2[1]) + ":" + str(freq_high2[1]) + "\t")
        outfile.write("frequency_fasta" + fname + "_" + str(unit2[1]) + ":" + str(freq_low2[1]) + "\t")
        outfile.write("relative_frequency")
        for i in range(len(iunit2)):
            outfile.write("\tfrequency_fasta" + fname + "_" + str(iunit2[i]) + "_in_"  + str(unit) + "_" + str(freq_high))
            outfile.write("\tfrequency_fasta" + fname + "_" + str(iunit2[i]) + "_in_" + str(unit) + "_" + str(freq_low))
            outfile.write("\trelative_frequency")
        outfile.write("\n")
        for i in range(len(list_high_freq[0])):
            for j in range(len(list_high_freq)):
                outfile.write(str(list_high_freq[j][i]) + "\t" + str(list_low_freq[j][i]) + "\t" + str(relative_freq[j][i]) + "\t")
            outfile.write("\n")
        for i in range(len(mean_high)):
            outfile.write(str(mean_high[i]) + "\t" + str(mean_low[i]) + "\t" + str(mean_rel_freq[i]) + "\t")
        outfile.write("mean\n")
        for i in range(len(std_high)):
            outfile.write(str(std_high[i]) + "\t" + str(std_low[i]) + "\t" + str(std_rel_freq[i]) + "\t")
        outfile.write("std\n")
        for i in range(len(pval)):
            outfile.write(str(pval[i]) + "\t \t \t ")
        outfile.write("pvalue\n")


def write_tsv_file(unit_type, unit, freq_high, freq_low, list_high_freq, list_low_freq, iteration, output):
    """
    Create a tsv file.

    :param unit_type: (string) the unit type we want to study : it can be *dnt* or *feature or *nt*
    :param unit: (string) the unit for which we want to calculate the median frequency
    :param freq_high: (float) the frequency of ``unit`` for the fasta file having an high content of the unit of \
    interest
        :param freq_low: (float) the frequency of ``unit`` for the fasta file having a low content of the unit of \
    interest
    :param list_high_freq: (list of float) list of mean frequencies in a particular unit in a list of fasta file enriched \
    for this unit
    :param list_low_freq: (list of float) list of mean frequencies in a particular unit in a list of fasta file impoverished \
    for this unit
    :param iteration: (int) the number of fasta files we want to create
    :param output: (string) path where the fasta_file will be created
    """
    unit2 = unit.split(",")
    fname = "mutated"
    freq_high2 = freq_high.split(",")
    freq_low2 = freq_low.split(",")
    relative_freq = [get_relative_freq_values(list_high_freq[0], list_low_freq[0]), get_relative_freq_values(list_high_freq[1], list_low_freq[1])]
    mean_rel_freq = [np.mean(relative_freq[0]), np.mean(relative_freq[1])]
    sd_rel_freq = [np.std(relative_freq[0], ddof=1), np.std(relative_freq[1], ddof=1)]

    mean_high = [np.mean(list_high_freq[0]), np.mean(list_high_freq[1])]
    std_high = [np.std(list_high_freq[0], ddof=1), np.std(list_high_freq[1], ddof=1)]

    mean_low = [np.mean(list_low_freq[0]), np.mean(list_low_freq[1])]
    std_low = [np.std(list_low_freq[0], ddof=1), np.std(list_low_freq[1], ddof=1)]
    unit = unit.replace("_", "-")
    with open(output + unit_type + "_" + unit + "_frequency_comparison_between_" + str(iteration) + "_" + fname + "_fasta_file-high:"+str(freq_high) + "_low:" + str(freq_low) + ".tsv", "w") as outfile:
        outfile.write("frequency_fasta" + fname + "_" + str(unit2[0]) + ":" + str(freq_high2[0]) + "\t")
        outfile.write("frequency_fasta" + fname + "_" + str(unit2[0]) + ":" + str(freq_low2[0]) + "\t")
        outfile.write("relative_frequency\t")
        outfile.write("frequency_fasta" + fname + "_" + str(unit2[1]) + ":" + str(freq_high2[1]) + "\t")
        outfile.write("frequency_fasta" + fname + "_" + str(unit2[1]) + ":" + str(freq_low2[1]) + "\t")
        outfile.write("relative_frequency\n")
        for i in range(len(list_high_freq[0])):
            outfile.write(str(list_high_freq[0][i]) + "\t" + str(list_low_freq[0][i]) + "\t" + str(relative_freq[0][i]) + "\t")
            outfile.write(str(list_high_freq[1][i]) + "\t" + str(list_low_freq[1][i]) + "\t" + str(relative_freq[1][i]) + "\n")
        outfile.write(str(mean_high[0]) + "\t" + str(mean_low[0]) + "\t" + str(mean_rel_freq[0]) + "\t" +
                      str(mean_high[1]) + "\t" + str(mean_low[1]) + "\t" + str(mean_rel_freq[1]) + "\t" + "mean\n")
        outfile.write(str(std_high[0]) + "\t" + str(std_low[0]) + "\t" + str(sd_rel_freq[0]) + "\t" +
                      str(std_high[1]) + "\t" + str(std_low[1]) + "\t" + str(sd_rel_freq[1]) + "\t" +"std")


def fusion_stretch_count_list(stretch_count_up, stretch_count_down, stretch_count_rel):
    """
    Fusion each sublist in stretch_count_up and stretch_count_down one after the other!

    :param stretch_count_up: (list of list of int) give the proportion of exons having at least one strecth \
    of size given by global variable ``strecthes`` (for each stretch in strecthes there's a sublist in \
    stretch_count_up)
    :param stretch_count_down: (list of list of int) give the proportion of exons having at least one strecth \
    of size given by global variable ``strecthes`` (for each stretch in strecthes there's a sublist in \
    stretch_count_down)
    :param stretch_count_rel: (list of list of float) give the relative proportion of exons having at least 1 \
    streact of the different kind of stretch studied (each sublist is the relative frequency of a \
    stretch for the given number of fasta file created)
    :return: A list containing a sublist of stretch_count_up followed by a sublist of stretch_count_down
    """
    merged_list = []
    name_list = []
    for i in range(len(stretch_count_up)):
        cur_stretch = str(stretches[i][0]) + "/" + str(stretches[i][1])
        merged_list.append(stretch_count_up[i])
        name_list.append(cur_stretch + "_high")
        merged_list.append(stretch_count_down[i])
        name_list.append(cur_stretch + "_low")
        merged_list.append(stretch_count_rel[i])
        name_list.append(cur_stretch + "_relative_freq")
    return merged_list, name_list


def write_stretch_table(strecth_count_list, name_strecth_list, size_fasta, filename_high, filename_low, output, iunit_type, iunit, iteration):
    """

    :param strecth_count_list: (list of list of int) each sublist correspond to a fasta file, each number \
    in the sublist is the number of random sequences in the fasta file having at least one strecth of the \
    different lind of strecthes of interest define in global variable ``stretches``
    :param name_stretch_list: (string) link each streches name to their value in strecth_count_list
    :param size_fasta: (int) the number of sequence in every fasta file created
    :param filename_high: (string) the name of the fasta file based on the chosen features and its high proportion used  \
    to created it
    :param filename_low: (string) the name of the fasta file based on the chosen features and its low proportion used  \
    to created it
    :param output: (string) path where the file will be created
    """
    outfile = "STRETCH_of_" + str(iunit_type) + "_" + iunit + "_in_" + str(iteration) + "_fasta_high;" + \
              str(filename_high) + "_fasta_low:" + filename_low + ".tsv"

    means = [np.mean(strecth_count_list[i]) * 100 for i in range(len(strecth_count_list))]
    stds = [np.std(strecth_count_list[i], ddof=1) * 100 for i in range(len(strecth_count_list))]
    pvalues = ""
    for i in range(len(strecth_count_list)):
        if i%3 == 0:
            pvalues += str(stretch_evalutator.r_ttest(strecth_count_list[i], strecth_count_list[i+1])) + "\t"
        else:
            pvalues += " \t"
    pvalues += "t-test pvalue"


    with open(output + outfile, "w") as myfile:
        myfile.write("This table contains the proportion of sequences having at least one stretch of " + iunit_type + " " + iunit + "\n")
        myfile.write("The sequences was generated in two type of fasta files : "+ "\n")
        myfile.write(str(iteration) + " 'high' fasta file  : " + str(filename_high) + "\n")
        myfile.write(str(iteration) +  " 'low' fasta file  : " + str(filename_low) + "\n")
        myfile.write("Every fasta files contains " + str(size_fasta[0]) + " sequences" "\n")
        header = "\t".join(name_strecth_list)
        myfile.write(header + "\n")
        for i in range(len(strecth_count_list[0])):
            line = ""
            for j in range(len(strecth_count_list)):
                line += str(round(strecth_count_list[j][i] * 100,4)) + "\t"
            myfile.write(line[:-2] + "\n")

        myfile.write("\t".join(map(str, means)) + "\tMean\n")
        myfile.write("\t".join(map(str, stds)) + "\tStd\n")
        myfile.write(pvalues)

def relative_frequency_stretch(stretch_count_high, stretch_count_low):
    """

    :param stretch_count_high: (list of list of float) each sublist correspond to the proportion of exons having a particular type of \
    stretch in a fasta
    :param stretch_count_low: (list of list of float) each sublist correspond to the proportion of exons having a particular type of \
    stretch in a fasta
    :return:
        * The relative frequency of of stretch between the count up and the count down for the different kind of \
        stretch studied

    """
    rel_freq = [[] for j in range(len(stretches))]
    frac_rel = [[] for j in range(len(stretches))]
    for j in range(len(stretches)):
        for i in range(len(stretch_count_high[j])):
            try:
                rel_freq[j].append((float(stretch_count_high[j][i] - stretch_count_low[j][i])) / stretch_count_low[j][i])
                frac_rel[j].append(str(round(rel_freq[j][-1],4)))
            except ZeroDivisionError:
                rel_freq[j].append(float("nan"))
                frac_rel[j].append("nan")
    return rel_freq, frac_rel



def main(type_unit, unit, freq_high, freq_low, iteration, output, iscub, itype_unit, iunit):
    """
    Create file file that the frequencies of ``unit`` of many (``ieration``) random fasta enriched for this unit to \
    the frequencies for that ``unit`` of many (``ieration``) random fasta imoverished for this unit

    :param type_unit: (string) the unit type we want to enriched/impoverish : it can be *dnt* or *feature or *nt*
    :param unit: (string) the unit for which we want to calculate the median frequency
    :param freq_high: (float) the frequency of ``unit`` for the fasta file having an high content of the unit of \
    interest
        :param freq_low: (float) the frequency of ``unit`` for the fasta file having a low content of the unit of \
    interest
    :param iteration: (int) the number of fasta files we want to create
    :param output: (string) path where the fasta_file will be created
    :param iscub: (bollean) True if we want to create random fasta respecting CCE codon bias usage \
    False if we want to mutate fasta sequence
    :param itype_unit:  (string) the unit type we want to study, it can be *dnt* or *feature or *nt*
    :param iunit: the unit of ``itype_unit`` for which we want to calculate the mean frequency in fasta files
    """
    if itype_unit is None:
        list_high_freq = get_mean_frequency_of_multiple_fasta(type_unit, unit, freq_high, iteration, output,
                                                              itype_unit, iunit)
        list_low_freq = get_mean_frequency_of_multiple_fasta(type_unit, unit, freq_low, iteration, output,
                                                             itype_unit, iunit)
        write_tsv_file(type_unit, unit, freq_high, freq_low, list_high_freq, list_low_freq, iteration, output)
    else:
        list_high_freq = get_mean_frequency_of_multiple_fasta(type_unit, unit, freq_high, iteration, output,
                                                              itype_unit, iunit)
        list_low_freq = get_mean_frequency_of_multiple_fasta(type_unit, unit, freq_low, iteration, output,
                                                             itype_unit, iunit)
        if itype_unit == "nt" or itype_unit.split(",")[0] == "nt":
            strech_count_high= list_high_freq[-2]
            strech_count_low = list_low_freq[-2]
            strech_rel_val, strech_rel_frac = relative_frequency_stretch(strech_count_high, strech_count_low)
            size_low = list_low_freq[-1]
            size_high = list_high_freq[-1]
            if size_low != size_high:
                print("Warning : the low and the high fasta did not have the same number of sequences")
        list_high_freq = list_high_freq[:-2]
        list_low_freq = list_low_freq[:-2]
        write_full_tsv_file(type_unit, unit, freq_high, freq_low, list_high_freq, list_low_freq, iteration, output, itype_unit, iunit)

        if itype_unit == "nt" or itype_unit.split(",")[0] == "nt":
            itype_unit = itype_unit.split(",")[0]
            iunit = iunit.split(",")[0]
            merged_list, name_list = fusion_stretch_count_list(strech_count_high, strech_count_low, strech_rel_val)
            filename_high = str(unit) + "_" + str(freq_high)
            filename_low = str(unit) + "_" + str(freq_low)
            write_stretch_table(merged_list, name_list, size_low, filename_high, filename_low, output,
                                itype_unit, iunit, iteration)
            outfile = "STRETCH_of_" + str(itype_unit) + "_" + iunit + "_in_" + str(iteration) + "_fasta_high;" + \
                      str(filename_high) + "_fasta_low:" + filename_low + ".png"

            name_file = ["high : " + str(filename_high), "low : " + str(filename_low)]

            stretch_name = []
            for stretch_content, stretch_len in stretches:
                stretch_name.append(str(stretch_content) + "/" + str(stretch_len))

            stretch_evalutator.barplot_maker(merged_list, stretch_name, name_file, itype_unit, iunit, output + outfile, size_low[0],
                                             iteration)



def launcher():
    """
    function that contains a parser to launch the program
    """
    # description on how to use the program
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description="""
The goal of this script is to:
     1. Create a lot of random fasta enriched in a particular unit (*nucleotide*, *di-nucleotide* or *feature*)
     2. Create a lot of random fasta impoverished in **the same unit** (*nucleotide*, *di-nucleotide* or *feature*)

And then compare their frequency for this unit for each couple of random fasta enriched and impoverished \
(for this unit). To compare their frequency, the relative frequency is computed

The relative frequency is calculated as follow:

  F_{relative} = frac{F_{interest} - F_{control}}{F_{control}}

Where:
  * :math:`F_{relative}` is the relative frequency of a unit :math:`F`
  * :math:`F_{interest}` is the frequency of :math:`F` in the interest set of exons
  * :math:`F_{control}` is the frequency of :math:`F` in the control sets of exons

We then calculate the mean and the standard deviation of:
    1. The list of frequencies in **unit** obtained by creating random fasta files enriched in this **unit**
    2. The list of frequencies in **unit** obtained by creating random fasta files impoverished in this **unit**
    3. The list of relative frequencies between random fasta files enriched and impoverished in this **unit**
    """)
    # Arguments for the parser

    parser.add_argument('--iteration', dest='iteration', help="the number of random sequence to cerate",
                        default=100)

    parser.add_argument('--iscub', dest='iscub', help="True to create cub sequence False to create mutated sequence",
                        default=False)
    parser.add_argument('--type_unit_interest', dest='type_unit_interest', help="unit type for which we want to calculate the average frequency in fasta files",
                        default=None)
    parser.add_argument('--unit_interest', dest='unit_interest', help="unit for which we want to calculate the average frequency in fasta files",
                        default=None)
    req_arg = parser.add_argument_group("required arguments")

    req_arg.add_argument('--type_unit', dest='type_unit', help="the type of unit we want to enrich (nt, dnt or feature)",
                         required=True)
    req_arg.add_argument('--unit', dest='unit', help="the unit for which we want to create random fasta having a low "
                                                     "content of if and a high content of it",
                         required=True)
    req_arg.add_argument('--freq_high', dest='freq_high', help="the frequency of unit in your enriched fasta file "
                                                              "enriched for this unit",
                        required=True)
    req_arg.add_argument('--freq_low', dest='freq_low', help="the frequency of unit in your impoverished fasta file "
                                                              "enriched for this unit",
                        required=True)
    req_arg.add_argument('--output', dest='output', help="path where the result will be created",
                        required=True)



    args = parser.parse_args()

    if args.output[-1] != "/":
        args.output += "/"
    if not os.path.isdir(args.output):
        print("ERROR : the given output folder does not exist")
        exit(1)

    if args.iscub == "False":
        args.iscub = False

    if args.iscub == "True":
        args.iscub = True


    if args.type_unit not in ["nt", "dnt", "feature"] and args.type_unit_interest not in ["nt", "dnt", "feature", None]:
        print("Unrecognized type unit")
        exit(1)

    try:
        args.iteration = int(args.iteration)
    except ValueError:
        print("Wrong value for iteration argument...")
        exit(1)

    main(args.type_unit, args.unit, args.freq_high, args.freq_low, args.iteration, args.output, args.iscub, args.type_unit_interest, args.unit_interest)

if __name__ == "__main__":
    launcher()
