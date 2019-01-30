# coding=utf-8
#################################################################
#                      Description
#################################################################
# script to make all the figures of the tRNA

#################################################################
#                        Imports
#################################################################


from dictionnary import codon2rareness
from dictionnary import codon2aminoAcid
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
from scipy.stats import *
import numpy as np
import sys
file_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, file_dir + "/frequence_recap_graphics/")


################################################################
#                      Function
################################################################


def create_an_exon_type_frequency_dictionary(exon_type):
    """
    :param exon_type: (string) the type of the exon of the control sets
    :return: 2 dictionaries : one containing the frequencies of every codon in fasterDB the other containing
    the frequencies of every amino acid in fasterDB
    """
    # saving the path
    file_dir = os.path.dirname(os.path.realpath(__file__))
    path_file_aa = file_dir + "/frequences/frequence/aa_frequency_" \
                   + str(exon_type) + ".txt"
    path_file_codon = file_dir + "/frequences/frequence/codon_frequency_" \
                      + str(exon_type) + ".txt"
    # creating the dictionaries to fill
    aa_dic = {}
    codon_dic = {}
    # filling the aa_dic
    amino_acid_file = open(path_file_aa, "r")
    amino_acid_file.readline()
    line = amino_acid_file.readline()
    while line:
        line = line.split("\t")
        aa_dic[line[0]] = float(line[1])
        line = amino_acid_file.readline()
    amino_acid_file.close()
    # filling the codon_dic
    codon_file = open(path_file_codon, "r")
    codon_file.readline()
    line = codon_file.readline()
    while line:
        line = line.split("\t")
        codon_dic[line[0]] = float(line[1])
        line = codon_file.readline()
    codon_file.close()
    return codon_dic, aa_dic


def create_a_corrected_pvalue_dictionary(output, type_unit):
    """
    :param output: (string), the path were to read the files that contains the corrected p-value
    :param type_unit: (string) the type of the unit for which we want to create a dictionary containing the
    corrected value of each unit of this type
    :return: a dictionary with the corrected p-value
    """
    p_corrected = {}
    if type_unit == "codon":
        file_name = "codon_enrichment_file.csv"
    else:
        file_name = "amino_acid_enrichment_file.csv"
    my_file = open(output + file_name, "r")
    my_file.readline()
    line = my_file.readline()
    while line:
        line = line.split(';')
        if type_unit == "codon":
            p_corrected[line[0]] = float(line[7])
        else:
            p_corrected[line[0]] = float(line[5])
        line = my_file.readline()
    return p_corrected


def cm2inch(value):
    """
    :param value: (float) value in centimeter
    :return: (float) value in inch
    """
    return value / 2.54


def create_a_graphic_with_control(prop_dic, pvalue_dic, corrected_pvalue_dic, control_prop_dic, type_unit, exon_type,
                                  output, name):
    """
    :param prop_dic: (dictionary) the proportion of all codon/amino acid in the interest set of exons
    :param pvalue_dic: (dictionary) the p-value for each codon/amino acid of the interest set
    :param corrected_pvalue_dic: (dictionary) the fdr corrected p-value for each codon/amino acid of the interest set
    :param control_prop_dic: (dictionary) the proportion for each amino acid/codon of a certain type of all exon with
    the type exon_type (ACE/FCE/LCE/CCE/LCE) in fasterDB
    :param type_unit: (string) codons or acides amines
    :param exon_type: (string) the type of the exon of the control sets
    :param output: (string) the path were the graph will be created
    :param name: string corresponding to the name of the set of exon that will figure in the graphics produced by the
    program
    Create a plot of the frequencies in codon/amino_acid in the exon of the user set in function of the frequencies
    """

    fig, ax = plt.subplots(figsize=(cm2inch(48), cm2inch(27)))
    abscissa = list()
    ordinate = list()
    abscissa_red = list()
    ordinate_red = list()
    abscissa_pink = list()
    ordinate_pink = list()
    max_o = 0
    max_a = 0
    for key in prop_dic.keys():
        if prop_dic[key] > max_o:
            max_o = prop_dic[key]
        if control_prop_dic[key] > max_a:
            max_a = control_prop_dic[key]
    if type_unit == "codons":
        val = max_a / 75
    else:
        val = 0
    if max_o > max_a:
        ax.set(xlim=(0, max_o + 0.002), ylim=(0, max_o + 0.002))
    else:
        ax.set(xlim=(0, max_a + 0.002), ylim=(0, max_a + 0.002))
    ax.plot()
    for key in prop_dic.keys():
        if type_unit == "codons":
            ax.annotate(str(key) + "$^{" + codon2aminoAcid[key] + codon2rareness[key] + "}$",
                        xy=(control_prop_dic[key], prop_dic[key]),
                        xytext=(control_prop_dic[key] - val, prop_dic[key] + max_o / 65))
        else:
            ax.annotate(key, xy=(control_prop_dic[key], prop_dic[key]),
                        xytext=(control_prop_dic[key] - val, prop_dic[key] + max_o / 65))
        if corrected_pvalue_dic[key] <= 0.05:
            abscissa_red.append(control_prop_dic[key])
            ordinate_red.append(prop_dic[key])
        elif pvalue_dic[key] <= 0.05:
            abscissa_pink.append(control_prop_dic[key])
            ordinate_pink.append(prop_dic[key])
        else:
            abscissa.append(control_prop_dic[key])
            ordinate.append(prop_dic[key])
    ax.plot(abscissa_red, ordinate_red, 'ro', label="fdr<=0.05")
    ax.plot(abscissa_pink, ordinate_pink, color="pink", marker="o", linewidth=0, label="p<=0.05")
    ax.plot(abscissa, ordinate, 'ko', label="non significatifs")

    ax.plot(ax.get_xlim(), ax.get_ylim(), color="red")

    if exon_type == "ALL":
        exon_type = ""

    xlabel = u"Fréquences en " + unicode(type_unit) + u" dans tous les exons " + unicode(exon_type) + u" de fasterDB"
    ylabel = u"Fréquences en " + unicode(type_unit) + u" dans les exons " + unicode(name)
    title = u"Fréquences en " + unicode(type_unit) + u" dans les exons " + unicode(name) + \
            u" en fonction des fréquences en " + unicode(type_unit) + u" dans tous les exons " + unicode(exon_type) + \
            u" de fasterDB"
    ax.legend(loc="upper left", shadow=True, numpoints=1)
    plt.xlabel(xlabel.replace("  ", " "))
    plt.title(title.replace("  ", " "))
    plt.ylabel(ylabel.replace("  ", " "))
    plt.savefig(output + str(type_unit) + "_figure.pdf", bbox_inches='tight')
    plt.savefig(output + str(type_unit) + "_figure.png", bbox_inches='tight')
    plt.clf()
    plt.cla()


def create_graphs(codon_frequencies, aa_frequencies, dic_p_val_codon, dic_p_val_aa, p_adjust_codon, p_adjust_aa,
                  exon_type, output, name):
    """
    :param codon_frequencies: (dictionary) a dictionary that links for each codon, it's frequency, in the exon set
    given by the user
    :param aa_frequencies: (dictionary) a dictionary that links for each amino acid, it's frequency, in the exon set
    given by the user
    :param dic_p_val_codon: (dictionary) a dictionary that give for each codon its p-value (impoverishment/enrichment
    test)
    :param dic_p_val_aa: (dictionary) a dictionary that give for each amino acid its p-value (impoverishment/enrichment
    test)
    :param p_adjust_codon: (dictionary) the fdr corrected p-value for each codon acid of the interest set
    :param p_adjust_aa: (dictionary) the fdr corrected p-value for each amino acid acid of the interest set
    :param exon_type: (string) the type of exon that will be compared to the exon set given by the user
    :param output: (string) the path were the results of the program will be stored
    :param name: (string) the name of the set of exon that will figure in the graphics produced by the program
    Create 2 graphics, one indicating the codon frequencies in the exons of the interest set in function of the codon
    frequencies in all the exon of type "exon type" of fasterDB
    """
    name = name.strip()
    control_dic_codon, control_dic_aa = create_an_exon_type_frequency_dictionary(exon_type)
    os.mkdir(output + "graphics")
    full_output = output + "graphics/"
    create_a_graphic_with_control(codon_frequencies, dic_p_val_codon, p_adjust_codon, control_dic_codon, "codons",
                                  exon_type, full_output, name)
    create_a_graphic_with_control(aa_frequencies, dic_p_val_aa, p_adjust_aa, control_dic_aa, "acides amines",
                                  exon_type, full_output, name)


def create_a_graphic_up_vs_down(prop_dic_up, prop_dic_down, type_unit, output, name, dic_nt_up=None, dic_nt_down=None):
    """
    :param prop_dic_up: (dictionary) the proportion of all codon/amino acid in the interest set of exons up
    :param prop_dic_down: (dictionary) the proportion of all codon/amino acid in the interest set of exons down
    :param type_unit: (string) codons or acides amines
    :param output: (string), the path were the graphic will be created
    :param name: (string), the named of the set of exon used
    :param dic_nt_up: (dictionary of float) it contains for each nucleotide its proportion in the set of exon
    up-regulated give by the user
    :param dic_nt_down: (dictionary of float) it contains for each nucleotide its proportion in the set of exon
    down-regulated give by the user
    """
    name = name.strip()
    abscissa = list()
    ordinate = list()
    fig, ax = plt.subplots(figsize=(cm2inch(48), cm2inch(27)))
    max_o = 0
    max_a = 0
    for key in prop_dic_up.keys():
        if prop_dic_up[key] > max_o:
            max_o = prop_dic_up[key]
        if prop_dic_down[key] > max_a:
            max_a = prop_dic_down[key]
    if type_unit == "codons":
        val = max_a / 75
    else:
        val = 0
    if max_o > max_a:
        ax.set(xlim=(0, max_o + 0.002), ylim=(0, max_o + 0.002))
    else:
        ax.set(xlim=(0, max_a + 0.002), ylim=(0, max_a + 0.002))
    for key in prop_dic_up.keys():
        if type_unit == "codons":
            ax.annotate(key + "$^{" + codon2aminoAcid[key] + codon2rareness[key] + "}$",
                        xy=(prop_dic_down[key], prop_dic_up[key]),
                        xytext=(prop_dic_down[key] - val, prop_dic_up[key] + max_o / 65))
        else:
            ax.annotate(key, xy=(prop_dic_down[key], prop_dic_up[key]),
                        xytext=(prop_dic_down[key] - val, prop_dic_up[key] + max_o / 65))
        abscissa.append(prop_dic_down[key])
        ordinate.append(prop_dic_up[key])
    ax.plot(abscissa, ordinate, 'ko')
    ax.plot(ax.get_xlim(), ax.get_ylim(), color="red")

    if type_unit == "codons":
        nt_up = """
                A = """ + str(dic_nt_up['A']) + """ %
                T = """ + str(dic_nt_up['T']) + """ %
                C = """ + str(dic_nt_up['C']) + """ %
                G = """ + str(dic_nt_up['G']) + """ %
                """
        nt_down = """
                  A = """ + str(dic_nt_down['A']) + """ %
                  T = """ + str(dic_nt_down['T']) + """ %
                  C = """ + str(dic_nt_down['C']) + """ %
                  G = """ + str(dic_nt_down['G']) + """ %
                  """
        ax.text(0.1, 0.9, nt_up, ha='center', va='center', transform=ax.transAxes)
        ax.text(0.9, 0.1, nt_down, ha='center', va='center', transform=ax.transAxes)
    xlabel = u"Fréquences en " + unicode(type_unit) + u" dans les exons down " + unicode(name)
    ylabel = u"Fréquences en " + unicode(type_unit) + u" dans les exons up " + unicode(name)
    title = u"Fréquences en " + unicode(type_unit) + u" dans les exons up en fonction des fréquences en " + \
            unicode(type_unit) + u" dans tous exons down"
    plt.xlabel(xlabel.replace("  ", " "))
    plt.title(title.replace("  ", " "))
    plt.ylabel(ylabel.replace("  ", " "))

    plt.savefig(output + str(type_unit) + "_figure_up_vs_down.pdf", bbox_inches='tight')
    plt.savefig(output + str(type_unit) + "_figure_up_vs_down.png", bbox_inches='tight')

    plt.clf()
    plt.cla()


def graphic_up_vs_down_stat(prop_dic_up, prop_dic_down, type_unit, dic_p_value, dic_p_value_corrected, output, name):
    """
    :param prop_dic_up: (dictionary) the proportion of all codon/amino acid in the interest set of exons up
    :param prop_dic_down: (dictionary) the proportion of all codon/amino acid in the interest set of exons down
    :param type_unit: (string) codons or acides amines
    :param dic_p_value: (dictionary of float) dictionary of p-value for each features (key of the dictionary) - feature
    can be codon or amino acids
    :param dic_p_value: (dictionary of float) dictionary of corrected p-value for each features (key of the dictionary)
    - feature can be codon or amino acids
    :param output: (string), the path were the graphic will be created
    :param name: (string), the named of the set of exon used
    """
    name = name.strip()
    abscissa = list()
    ordinate = list()
    abscissa_red = list()
    ordinate_red = list()
    abscissa_pink = list()
    ordinate_pink = list()
    abscissa_yellow = list()
    ordinate_yellow = list()
    fig, ax = plt.subplots(figsize=(cm2inch(48), cm2inch(27)))
    max_o = 0
    max_a = 0
    for key in prop_dic_up.keys():
        cur_mean = np.mean(prop_dic_up[key])
        if cur_mean > max_o:
            max_o = cur_mean
        cur_mean = np.mean(prop_dic_down[key])
        if cur_mean > max_a:
            max_a = cur_mean
    if type_unit == "codons":
        val = max_a / 75
    else:
        val = 0
    if max_o > max_a:
        ax.set(xlim=(0, max_o + 0.002), ylim=(0, max_o + 0.002))
    else:
        ax.set(xlim=(0, max_a + 0.002), ylim=(0, max_a + 0.002))
    for key in prop_dic_up.keys():
        cur_mean_up = np.mean(prop_dic_up[key])
        cur_mean_down = np.mean(prop_dic_down[key])
        if type_unit == "codons":
            ax.annotate(key + "$^{" + codon2aminoAcid[key] + codon2rareness[key] + "}$",
                        xy=(cur_mean_down, cur_mean_up),
                        xytext=(cur_mean_down - val, cur_mean_up + max_o / 65))
        else:
            ax.annotate(key, xy=(cur_mean_down, cur_mean_up),
                        xytext=(cur_mean_down - val, cur_mean_up + max_o / 65))
        if dic_p_value_corrected[key] <= 0.05:
            abscissa_red.append(cur_mean_down)
            ordinate_red.append(cur_mean_up)
        elif dic_p_value[key] <= 0.05:
            abscissa_pink.append(cur_mean_down)
            ordinate_pink.append(cur_mean_up)
        elif dic_p_value[key] == "NaN":
            abscissa_yellow.append(cur_mean_down)
            ordinate_yellow.append(cur_mean_up)
        else:
            abscissa.append(cur_mean_down)
            ordinate.append(cur_mean_up)
    ax.plot(abscissa_red, ordinate_red, 'ro', label="fdr<=0.05")
    ax.plot(abscissa_pink, ordinate_pink, color="pink", marker="o", linewidth=0, label="p<=0.05")
    ax.plot(abscissa_yellow, ordinate_yellow, color="yellow", marker="o", linewidth=0, label="NaN")
    ax.plot(abscissa, ordinate, 'ko', label="non significatifs")
    ax.plot(ax.get_xlim(), ax.get_ylim(), color="red")

    xlabel = u"Fréquences en " + unicode(type_unit) + u" dans les exons down " + unicode(name)
    ylabel = u"Fréquences en " + unicode(type_unit) + u" dans les exons up " + unicode(name)
    title = u"Fréquences en " + unicode(type_unit) + u" dans les exons up en fonction des fréquences en " + \
            unicode(type_unit) + u" dans tous exons down (exons > 9 codons)"

    ax.legend(loc="upper left", shadow=True, numpoints=1)
    plt.xlabel(xlabel.replace("  ", " "))
    plt.title(title.replace("  ", " "))
    plt.ylabel(ylabel.replace("  ", " "))

    plt.savefig(output + str(type_unit) + "_figure_up_vs_down_stat.pdf", bbox_inches='tight')
    plt.savefig(output + str(type_unit) + "_figure_up_vs_down_stat.png", bbox_inches='tight')

    plt.clf()
    plt.cla()


def graphic_up_vs_down_param(dic_param, type_unit, output):
    """
    :param dic_param: (dictionary of list of 6 float) - each key of the dictionary corresponds to a codon/aa -
    the list corresponds to the estimation of 3 parameters by gamlss for up and down : [sigma_up, mu_up, nu_up,
    sigma_down, mu_down, nu_down)
    :param type_unit: (string) the type of the feature we are analysing : codons or amino acids
    :param output: (string) path were the graphics will be created
    """
    fig = plt.figure(figsize=(21. /2.54, 44.5/2.54 ))
    # sigma
    ax = plt.subplot2grid((3, 1), (0, 0))
    abscissa = list()
    ordinate = list()
    max_o = 0
    max_a = 0
    for key in dic_param.keys():
        cur_mean = float(dic_param[key][0])
        if cur_mean > max_o and cur_mean != 'NaN' and key not in ["*", "TGA", "TAA", "TAG"]:
            max_o = cur_mean
        cur_mean = float(dic_param[key][3])
        if cur_mean > max_a and cur_mean != 'NaN' and key not in ["*", "TGA", "TAA", "TAG"]:
            max_a = cur_mean
    max_a = float(max_a)
    max_o = float(max_o)
    if type_unit == "codons":
        val = max_a / 75
    else:
        val = 0
    if max_o > max_a:
        ax.set(xlim=(0, max_o + 0.002), ylim=(0, max_o + 0.002))
    else:
        ax.set(xlim=(0, max_a + 0.002), ylim=(0, max_a + 0.002))
    for key in dic_param.keys():

        if type_unit == "codons":
            ax.annotate(key + "$^{" + codon2aminoAcid[key] + codon2rareness[key] + "}$",
                        xy=(float(dic_param[key][3]), float(dic_param[key][0])),
                        xytext=(float(dic_param[key][3]) - val, float(dic_param[key][0]) + max_o / 65))
        else:
            ax.annotate(key, xy=(float(dic_param[key][3]), float(dic_param[key][0])),
                        xytext=(float(dic_param[key][3]) - val, float(dic_param[key][0]) + max_o / 65))
        abscissa.append(float(dic_param[key][3]))
        ordinate.append(float(dic_param[key][0]))

    ax.plot(abscissa, ordinate, 'ko')
    ax.plot(ax.get_xlim(), ax.get_ylim(), color="red")

    ax.set_xlabel(u"Paramêtre sigma estimé par gamlss pour le jeu d'exons down ")
    ax.set_ylabel(u"Paramêtre sigma estimé par gamlss pour le jeu d'exons up ")
    ax.set_title(u"Comparaison du paramêtre sigma estimé par gamlss dans les exons up et down")
    fig.add_subplot(ax)

    ax = None
    ax = plt.subplot2grid((3, 1), (1, 0))
    abscissa = list()
    ordinate = list()
    max_o = 0
    max_a = 0
    for key in dic_param.keys():
        cur_mean = float(dic_param[key][1])
        if cur_mean > max_o and key not in ["*", "TGA", "TAA", "TAG"]:
            max_o = cur_mean
        cur_mean = float(dic_param[key][4])
        if cur_mean > max_a and key not in ["*", "TGA", "TAA", "TAG"]:
            max_a = cur_mean
    max_a = float(max_a)
    max_o = float(max_o)
    if type_unit == "codons":
        val = max_a / 75
    else:
        val = 0
    if max_o > max_a:
        ax.set(xlim=(0, max_o + 0.002), ylim=(0, max_o + 0.002))
    else:
        ax.set(xlim=(0, max_a + 0.002), ylim=(0, max_a + 0.002))
    for key in dic_param.keys():

        if type_unit == "codons":
            ax.annotate(key + "$^{" + codon2aminoAcid[key] + codon2rareness[key] + "}$",
                        xy=(float(dic_param[key][4]), float(dic_param[key][1])),
                        xytext=(float(dic_param[key][4]) - val, float(dic_param[key][1]) + max_o / 65))
        else:
            ax.annotate(key, xy=(float(dic_param[key][4]), float(dic_param[key][1])),
                        xytext=(float(dic_param[key][4]) - val, float(dic_param[key][1]) + max_o / 65))
        abscissa.append(float(dic_param[key][4]))
        ordinate.append(float(dic_param[key][1]))

    ax.plot(abscissa, ordinate, 'ko')
    ax.plot(ax.get_xlim(), ax.get_ylim(), color="red")

    ax.set_xlabel(u"Paramêtre mu estimé par gamlss pour le jeu d'exons down")
    ax.set_ylabel(u"Paramêtre mu estimé par gamlss pour le jeu d'exons up ")
    ax.set_title(u"Comparaison du paramêtre mu estimé par gamlss dans les exons up et down")
    fig.add_subplot(ax)

    ax = None
    ax = plt.subplot2grid((3, 1), (2, 0))
    abscissa = list()
    ordinate = list()
    max_o = 0
    max_a = 0
    for key in dic_param.keys():
        cur_mean = float(dic_param[key][2])
        if cur_mean > max_o and key not in ["*", "TGA", "TAA", "TAG"]:
            max_o = cur_mean
        cur_mean = float(dic_param[key][5])
        if cur_mean > max_a and key not in ["*", "TGA", "TAA", "TAG"]:
            max_a = cur_mean
    max_a = float(max_a)
    max_o = float(max_o)
    if type_unit == "codons":
        val = max_a / 75
    else:
        val = 0
    if max_o > max_a:
        ax.set(xlim=(0, max_o + 0.002), ylim=(0, max_o + 0.002))
    else:
        ax.set(xlim=(0, max_a + 0.002), ylim=(0, max_a + 0.002))
    for key in dic_param.keys():

        if type_unit == "codons":
            ax.annotate(key + "$^{" + codon2aminoAcid[key] + codon2rareness[key] + "}$",
                        xy=(float(dic_param[key][5]), float(dic_param[key][2])),
                        xytext=(float(dic_param[key][5]) - val, float(dic_param[key][2]) + max_o / 65))
        else:
            ax.annotate(key, xy=(float(dic_param[key][5]), float(dic_param[key][2])),
                        xytext=(float(dic_param[key][5]) - val, float(dic_param[key][2]) + max_o / 65))
        abscissa.append(float(dic_param[key][5]))
        ordinate.append(float(dic_param[key][2]))

    ax.plot(abscissa, ordinate, 'ko')
    ax.plot(ax.get_xlim(), ax.get_ylim(), color="red")

    ax.set_xlabel(u"Paramêtre nu estimé par gamlss pour le jeu d'exons down ")
    ax.set_ylabel(u"Paramêtre nu estimé par gamlss pour le jeu d'exons up ")
    ax.set_title(u"Comparaison du paramêtre nu estimé par gamlss dans les exons up et down")
    fig.add_subplot(ax)

    plt.savefig(output + str(type_unit) + "_up_vs_down_param.pdf")
    plt.savefig(output + str(type_unit) + "_up_vs_down_param.png")

    plt.clf()
    plt.cla()



def shap_var_maker(list_down, list_up):
    """
    :param list_down: (list of float) list of proportion of a nucleotide in every exon down-regulated given by the user
    :param list_up: (list of float) list of proportion of a nucleotide in every exon up-regulated given by the user
    :return: (boolean) true if the values in the list are normally distributed (p-value of the shapiro test < 0.05)
    and false else
    """
    x, p_val_down = shapiro(list_down)
    x, p_val_up = shapiro(list_up)

    if p_val_down > 0.05 and p_val_up > 0.05:
        fisher = np.var(list_up) / np.var(list_down)
        df_up = len(list_up) - 1
        df_down = len(list_down) - 1
        p_val = f.cdf(fisher, df_up, df_down)
        return p_val > 0.05
    else:
        return False


def make_label(label, list_up, list_down):
    """
    :param label: (string) a nucleotide letter (A,T,G,C) or 2 nucleotides letters or "purine" or "pyrimidine"
    :param list_up: (list of float) list of proportion of a nucleotide in every exon down-regulated given by the user
    :param list_down: (list of float) list of proportion of a nucleotide in every exon up-regulated given by the user
    :return: the p_value of the appropriate mean test comparison (equal variance and normality of the data for t-test
    else its a mann whitney test that is performed)
    """
    if shap_var_maker(list_down, list_up):
        test_stat = 't test'
        t_stat, p_val = ttest_ind(list_up, list_down)
    else:
        test_stat = 'mann whitney'
        u, p_val = mannwhitneyu(list_up, list_down)
    return label + " - " + "test " + test_stat + "  : p = " + str(p_val) + "\n"


def sum_list(nt_exons_up, nt_exons_down):
    """
    :param nt_exons_up: (list of floats) the frequencies of an nucleotide for each exons in the user's set of the
    up-reagulated exons
    :param nt_exons_down: (list of floats) the frequencies of an nucleotide for each exons in the user's set of the
    down-regulated exons
    :return:
    """
    res = []
    for i in range(len(nt_exons_up)):
        res.append(nt_exons_up[i] + nt_exons_down[i])
    return res


def create_an_up_vs_down_comparison(nt_exons_up, nt_exons_down, output):
    """
    :param nt_exons_up: (dictionary of list of floats) : For each nucleotides gives a list of the proportion of this
    nucleotide for each exons up-regulated given by the user
    :param nt_exons_down: (dictionary of list of floats) :  For each nucleotides gives a list of the proportion of this
    nucleotide for each exons down -regulated given by the user
    :param output: (string), the path where the graphic will be created
    """
    fig, ax = plt.subplots(figsize=(48. / 2.54, 27 / 2.54))
    labels = ['A\nup', 'A\ndown', 'T\nup', 'T\ndown', 'G\nup', 'G\ndown', 'C\nup', 'C\ndown',
              'Purine\nup', 'Purine\ndown', 'Pyrimidine\nup', 'Pyrimidine\ndown', 'GC\nup', 'GC\ndown',
              'AT\nup', 'AT\ndown']
    colors = sorted(['black', 'blue', 'orange', 'red', 'purple', 'green', 'brown', 'yellow'] * 2)
    box = ax.boxplot([nt_exons_up['A'], nt_exons_down['A'],  # Adenine comparison
                      nt_exons_up['T'], nt_exons_down['T'],  # Thymine comparison
                      nt_exons_up['G'], nt_exons_down['G'],  # Guanine comparison
                      nt_exons_up['C'], nt_exons_down['C'],  # Cytosine comparison
                      # Purine Comparison
                      sum_list(nt_exons_up['A'], nt_exons_up['G']), sum_list(nt_exons_down['A'], nt_exons_down['G']),
                      # Pyrimidine
                      sum_list(nt_exons_up['T'], nt_exons_up['C']), sum_list(nt_exons_down['T'], nt_exons_down['C']),
                      # GC Comparison
                      sum_list(nt_exons_up['G'], nt_exons_up['C']), sum_list(nt_exons_down['G'], nt_exons_down['C']),
                      # AT Comparison
                      sum_list(nt_exons_up['T'], nt_exons_up['A']), sum_list(nt_exons_down['T'], nt_exons_down['A']),
                      ], labels=labels, showmeans=True)

    val = 0
    for b in box['boxes']:
        b.set_color(colors[val])
        val += 1
    my_text = ""
    if len(nt_exons_up['A']) > 3:  # for shapiro to work with small list of exon
        my_text = "up vs down\n"
        my_text += make_label("A", nt_exons_up['A'], nt_exons_down['A'])
        my_text += make_label("T", nt_exons_up['T'], nt_exons_down['T'])
        my_text += make_label("G", nt_exons_up['G'], nt_exons_down['G'])
        my_text += make_label("C", nt_exons_up['C'], nt_exons_down['C'])
        my_text += make_label("Purine", nt_exons_up['A'] + nt_exons_up['G'], nt_exons_down['A'] + nt_exons_down['G'])
        my_text += make_label("Pyrimidine", nt_exons_up['T'] + nt_exons_up['C'], nt_exons_down['T'] +
                              nt_exons_down['C'])
        my_text += make_label("GC", nt_exons_up['G'] + nt_exons_up['C'], nt_exons_down['G'] + nt_exons_down['C'])
        my_text += make_label("AT", nt_exons_up['T'] + nt_exons_up['A'], nt_exons_down['T'] + nt_exons_down['A'])

    ax.text(0.9, 0.9, my_text, ha='center', va='center', transform=ax.transAxes)

    plt.title(u"Comparaison des proportions en nucléotides entre les ORF des exons up et down")
    plt.xlabel(u"Nucléotides / Groupes de nucléotides")
    plt.ylabel(u"Fréquences dans les exons")
    plt.savefig(output + "boxplot_nt_proportion_figure.pdf", bbox_inches='tight')
    plt.savefig(output + "boxplot_nt_proportion_figure.png", bbox_inches='tight')
    plt.clf()
    plt.cla()
    plt.close()


def remove_zero(prop_list):
    """
    :param prop_list: (list of float) list of proportion of a given codon, amino acid in a given set of exon
    :return: prop_list with all the zero proportion removed
    """
    new_list = []
    for value in prop_list:
        if value != 0:
            new_list.append(value)
    return new_list


def boxplot_composition(output, up_prop_list, down_prop_list, up_down_prop_list, unit_name):
    """
    Create a boxplot of the frequencies of a given amino acid/codon of the user up and down set of exon. It also
    create 2 other boxplot  of the frequencies of this amino acid/codon of the fasterDB ACE and CCE set of exon
    :param output:  (string) the output path where the figure will be created
    :param up_prop_list: (list of floats) list of proportion of a given codon, amino acid (i.e. information
    contained in the variable unit_name) in the up set of exon
    :param down_prop_list: (list of floats) list of proportion of a given codon, amino acid (i.e. information
    contained in the variable unit_name) in the down set of exon
    :param up_down_prop_list:(list of floats) list of proportion of a given codon, amino acid (i.e. information
    contained in the variable unit_name) in the merged up and down set of exons given by the user
    :param unit_name: (string) the name of the amino acid/ codon of interest
    """
    fig, ax = plt.subplots(figsize=(48. / 2.54, 27 / 2.54))
    # loading the variable we need to complete the graphics
    box_mod = __import__("result_ACE_CCE_boxplot_recap")

    if up_down_prop_list is not None:
        labels = ["$up \geq 0$", "$up > 0$", "$down \geq 0$", "$down > 0$", "$up\ down \geq 0$", "$up\ down > 0$",
                  "$ACE \geq 0$", "$ACE > 0$", "$CCE \geq 0$", "$CCE > 0$"]
        colors = sorted(['#0000FF', '#6C6CFF', '#FF6F00', '#F79F52', '#00FF09',
                         '#6EFFA8', "#EF00FF", "#F887FF", "#7A2900", "#C66B3E"])
    elif up_down_prop_list is None and down_prop_list is not None:
        labels = ["$up \geq 0$", "$up > 0$", "$down \geq 0$", "$down > 0$",
                  "$ACE \geq 0$", "$ACE > 0$", "$CCE \geq 0$", "$CCE > 0$"]
        colors = sorted(['#0000FF', '#6C6CFF', '#FF6F00', '#F79F52', "#EF00FF", "#F887FF", "#7A2900", "#C66B3E"])
    else:
        labels = ["$input \geq 0$", "$input > 0$", "$ACE \geq 0$", "$ACE > 0$", "$CCE \geq 0$", "$CCE > 0$"]
        colors = sorted(['#0000FF', '#6C6CFF', "#EF00FF", "#F887FF", "#7A2900", "#C66B3E"])

    # defining the variable that we will use
    ace_unit = [0]
    cce_unit = [0]

    # replacing stop aa names '*' by 'stop
    if unit_name == "*":
        unit_name = "stop"

    # To get the variable we want to use
    exec ("ace_unit =  box_mod.ACE_" + str(unit_name))
    exec ("cce_unit =  box_mod.CCE_" + str(unit_name))

    # removing zeros from the first list
    up_prop_list_without_zero = remove_zero(up_prop_list)
    down_prop_list_without_zero = None
    if down_prop_list is not None:
        down_prop_list_without_zero = remove_zero(down_prop_list)
    up_down_prop_list_without_zero = None
    if up_down_prop_list is not None:
        up_down_prop_list_without_zero = remove_zero(up_down_prop_list)
    ace_unit_without_zero = remove_zero(ace_unit)
    cce_unit_without_zero = remove_zero(cce_unit)

    if up_down_prop_list is not None:
        box = ax.boxplot([up_prop_list, up_prop_list_without_zero,
                          down_prop_list, down_prop_list_without_zero,
                          up_down_prop_list, up_down_prop_list_without_zero,
                          ace_unit, ace_unit_without_zero,
                          cce_unit, cce_unit_without_zero], labels=labels, showmeans=True)
    elif up_down_prop_list is None and down_prop_list is not None:
        box = ax.boxplot([up_prop_list, up_prop_list_without_zero,
                          down_prop_list, down_prop_list_without_zero,
                          ace_unit, ace_unit_without_zero,
                          cce_unit, cce_unit_without_zero], labels=labels, showmeans=True)
    else:
        box = ax.boxplot([up_prop_list, up_prop_list_without_zero,
                          ace_unit, ace_unit_without_zero,
                          cce_unit, cce_unit_without_zero], labels=labels, showmeans=True)

    val = 0
    for b in box['boxes']:
        b.set_color(colors[val])
        val += 1

    my_text = "information exon without zero\n"
    my_text += "up_list = " + str(len(up_prop_list) - up_prop_list.count(0)) + "/" + str(
        len(up_prop_list)) + " - " + str(
        round((1. - float(up_prop_list.count(0)) / len(up_prop_list)) * 100, 2)) + " %\n"
    if down_prop_list is not None:
        my_text += "down_list = " + str(len(down_prop_list) - down_prop_list.count(0)) + "/" + str(
                    len(down_prop_list)) + " - " + str(round((1. - float(down_prop_list.count(0)) /
                                                              len(down_prop_list)) * 100, 2)) + " %\n"
    if up_down_prop_list is not None:
        my_text += "up_down_list = " + str(len(up_down_prop_list) - up_down_prop_list.count(0)) + "/" + str(
                    len(up_down_prop_list)) + " - " + str(
                    round((1. - float(up_down_prop_list.count(0)) / len(up_down_prop_list)) * 100, 2)) + " %\n"
    my_text += "ACE_list = " + str(len(ace_unit) - ace_unit.count(0)) + "/" + str(
        len(ace_unit)) + " - " + str(
        round((1. - float(ace_unit.count(0)) / len(ace_unit)) * 100, 2)) + " %\n"
    my_text += "CCE_list = " + str(len(cce_unit) - cce_unit.count(0)) + "/" + str(
        len(cce_unit)) + " - " + str(
        round((1. - float(cce_unit.count(0)) / len(cce_unit)) * 100, 2)) + " %\n"

    ax.text(0.9, 0.9, my_text, ha='center', va='center', transform=ax.transAxes)
    """
    if down_prop_list is not None:
        my_text2 = "Test statistiques (sans correction) : \n"
        my_text2 += make_label("up vs down", up_prop_list, down_prop_list)
        my_text2 += make_label("up vs ACE", up_prop_list, ace_unit)
        my_text2 += make_label("up vs CCE", up_prop_list, cce_unit)
        my_text2 += make_label("down vs ACE", down_prop_list, ace_unit)
        my_text2 += make_label("down vs CCE", down_prop_list, cce_unit)
    else:
        my_text2 = "Test statistiques (sans correction) : \n"
        my_text2 += make_label("input vs ACE", up_prop_list, ace_unit)
        my_text2 += make_label("input vs CCE", up_prop_list, cce_unit)


    ax.text(0.25, 0.9, my_text2, ha='center', va='center', transform=ax.transAxes)
    """
    plt.title(u"Comparaison des proportions en " + unicode(unit_name) + u" entre les exons up et down et up+down")
    plt.xlabel(u"Jeux d'exons")
    plt.ylabel(u"Fréquences dans les exons")
    plt.savefig(output + str(unit_name) + "_boxplot_proportion_figure.png", bbox_inches='tight')
    plt.clf()
    plt.cla()
    plt.close()


def barplot_count_proportion(output, prop_count_list_up, prop_count_list_down, unit_name):
    """
    Function that creates 4 barplots per figure to get the proportion of exons that have 0, 1, 2, 3, 4, 5 or more
    occurrence of a given amino acid/ codons in four list of exon : the up regulated exon list, the down regulated exon
    list, the ACE exon list and the CCE exon list
    :param output: (string) the path where the graphics will be created
    :param prop_count_list_up: (list of six floats) list that gave the proportion of exons that have zero, one, two,
    three, four or five occurrence of a given codon/amino acid in the up exon list
    :param prop_count_list_down:(list of six floats) list that gave the proportion of exons that have zero, one, two,
    three, four or five occurrence of a given codon/amino acid in the down exon list
    :param unit_name: (string) the name of the amino acid/ codon of interest
    """
    fig = plt.figure()
    ax = plt.subplot2grid((2, 2), (0, 0))
    labels = ["0", "1", "2", "3+"]
    ax.bar([0, 1, 2, 3], [prop_count_list_up[0], prop_count_list_up[1], prop_count_list_up[2],
                                prop_count_list_up[3]],
           width=0.8, color="b")
    ax.set_xticks([0.4, 1.4, 2.4, 3.4])
    ax.set_xticklabels(labels)
    ax.set_title(u"Codon up : Proportion d'exons\navec $n$ " + unicode(unit_name) + u" $n \in [0;3]$")
    ax.grid()
    # ax.set_xlabel(u"Nombre de " + str(unit_name) + u" dans les exons" )
    ax.set_ylabel(u"Proportion d'exon ")
    fig.add_subplot(ax)

    if prop_count_list_down is not None:
        ax2 = plt.subplot2grid((2, 2), (0, 1))
        ax2.bar([0, 1, 2, 3], [prop_count_list_down[0], prop_count_list_down[1], prop_count_list_down[2],
                                     prop_count_list_down[3]],
                width=0.8, color="b", label=labels)
        ax2.set_xticks([0.4, 1.4, 2.4, 3.4])
        ax2.set_xticklabels(labels)
        ax2.set_title(u"Codon down : Proportion d'exons\navec $n$ " + unicode(unit_name) + u" $n \in [0;3]$")
        ax2.grid()
        fig.add_subplot(ax2)

    # ----------------------------------------- control figures ---------------------------------------------
    # loading the variable we need to complete the graphics
    bar_mod = __import__("result_ACE_CCE_barplot_recap")

    # defining the variable that we will use
    ace_unit = [0]
    cce_unit = [0]

    # replacing stop aa names '*' by 'stop
    if unit_name == "*":
        unit_name = "stop"

    # To get the variable we want to use
    exec("ace_unit =  bar_mod.ACE_" + str(unit_name))
    exec("cce_unit =  bar_mod.CCE_" + str(unit_name))

    # Creation of the ACE barplot
    ax3 = plt.subplot2grid((2, 2), (1, 0))
    print(ace_unit)
    ax3.bar([0, 1, 2, 3], [ace_unit[0], ace_unit[1], ace_unit[2], ace_unit[3]],
            width=0.8, color="b", label=labels)
    ax3.set_xticks([0.4, 1.4, 2.4, 3.4, 4.4, 5.4])
    ax3.set_xticklabels(labels)
    ax3.set_title(u"Proportion d'exons ACE\navec $n$ " + unicode(unit_name) + u" $n \in [0;5]$")
    ax3.set_xlabel(u"Nombre de " + str(unit_name) + u" dans les exons")
    ax3.set_ylabel(u"Proportion d'exon ")
    ax3.grid()
    fig.add_subplot(ax3)

    # Creation of the CCE barplot
    ax4 = plt.subplot2grid((2, 2), (1, 1))
    ax4.bar([0, 1, 2, 3], [cce_unit[0], cce_unit[1], cce_unit[2], cce_unit[3]],
            width=0.8, color="b", label=labels)
    ax4.set_xticks([0.4, 1.4, 2.4, 3.4, 4.4, 5.4])
    ax4.set_xticklabels(labels)
    ax4.set_title(u"Proportion d'exons CCE\navec $n$ " + unicode(unit_name) + u" $n \in [0;5]$")
    ax4.set_xlabel(u"Nombre de " + str(unit_name) + u" dans les exons")
    ax4.grid()
    # ax4.set_ylabel(u"Proportion d'exon ")
    fig.add_subplot(ax4)
    plt.tight_layout()
    plt.savefig(output + str(unit_name) + "_barplot_proportion_figure.pdf", bbox_inches='tight')
    plt.clf()
    plt.cla()
    plt.close()


def set_color(value, limit_value):
    """
    :param value: (int) the position of the boxplot of interest within the graph created by the functions
    boxplot_nucleotide_rich_or_poor_codon and boxplot_dinucleotide_rich_or_poor_codon
    :param limit_value: (int) the number of group of boxplot within the graph
    :return:
    """
    if value/limit_value == 0:
        return "blue"
    elif value/limit_value == 1:
        return "red"
    elif value/limit_value == 2:
        return "yellow"
    elif value / limit_value == 3:
        return "green"
    elif value / limit_value == 4:
        return "black"
    else:
        return "#00FFFF"


def boxplot_nucleotide_rich_or_poor_codon(output, up_dict, down_dict, up_down_dict, rich):
    """
    :param output: (string) the folder where the graphic will be created
    :param up_dict: (dictionary of list of float), a dictionary that gives for each up-regulated exons within a
    particular condition the proportion of codon rich/poor in (or that ends with)  A, C, G, or T nucleotides.
    A codon rich in A, for example, is a codon that contains at least 2 A.  A codon poor in A is a codon that contains
    no A
    :param down_dict: (dictionary of list of float or None if not present), a dictionary that gives for each
    down-regulated exons within a particular condition the proportion of codon rich/poor in (or that ends with) A, C, G,
     or T nucleotides. A codon rich in A, for example,is a codon that contains at least 2 A. A codon poor in A is a
     codon that contains no A
    :param up_down_dict: (dictionary of list of float or None if not present), a dictionary that gives for each
    regulated (up+down) exons within a particular condition the proportion of codon rich/poor in (or that ends with) A,
    C, G, or T nucleotides. A codon rich in A,  for example, is a codon that contains at least 2 A. A codon poor in A
    is a codon that contains no A
    :param rich: (string) "rich" if we want to create a graphics that gives the proportion of codons codons rich in A,
    C, G, T in the  up-down-up_and_down-regulated set of exons. "poor" if we want to create a graphics that gives the
    proportion of codons that contains no A, C, G, or T in the up-down-up_and_down-regulated set of exons. "last" if
    we want to create a graphics that gives the proportion of codons that ends with A, C, G, or T in the up-down-up_and_
    down-regulated set of exons
    """
    size_font = 20
    fig, ax = plt.subplots(figsize=(48. / 2.54, 27 / 2.54))
    # loading the variable we need to complete the graphics
    file_dir = os.path.dirname(os.path.realpath(__file__))
    sys.path.insert(0, file_dir + "/nucleotides_riche_codon_ACE_CCE/")
    if rich == "rich":
        A='A'; C='C'; G='G'; T='T'
        ACE_f = __import__("ACE_rich")
        CCE_f = __import__("CCE_rich")
    elif rich == "poor":
        A = 'B'; C = 'D'; G = 'H'; T = 'V'
        ACE_f = __import__("ACE_poor")
        CCE_f = __import__("CCE_poor")
    else:
        A = 'A'; C = 'C'; G = 'G'; T = 'T'
        ACE_f = __import__("ACE_last_nt")
        CCE_f = __import__("CCE_last_nt")
    dic_ACE = ACE_f.dic_ACE
    dic_CCE = CCE_f.dic_CCE
    if down_dict is not None and up_down_dict is not None:
        labels = ["up", "down", "up\nand\ndown", "ACE", "CCE"] * 4
        pos = [1, 2, 3, 4, 5, 7, 8, 9, 10, 11, 13, 14, 15, 16, 17, 19, 20, 21, 22, 23]
        ax.set_xticks(pos)
        box = ax.boxplot([up_dict[A], down_dict[A], up_down_dict[A], dic_ACE[A], dic_CCE[A],
                         up_dict[C], down_dict[C], up_down_dict[C], dic_ACE[C], dic_CCE[C],
                         up_dict[G], down_dict[G], up_down_dict[G], dic_ACE[G], dic_CCE[G],
                          up_dict[T], down_dict[T], up_down_dict[T], dic_ACE[T], dic_CCE[T]],
                         positions=pos, labels=labels, showmeans=True, notch=False, patch_artist=True)
        for i in range(len(box['boxes'])):
            box['boxes'][i].set_facecolor(set_color(i, 5))
        violin = ax.violinplot([up_dict[A], down_dict[A], up_down_dict[A], dic_ACE[A], dic_CCE[A],
                                up_dict[C], down_dict[C], up_down_dict[C], dic_ACE[C], dic_CCE[C],
                                up_dict[G], down_dict[G], up_down_dict[G], dic_ACE[G], dic_CCE[G],
                                up_dict[T], down_dict[T], up_down_dict[T], dic_ACE[T], dic_CCE[T]],
                               pos, showmeans=True, showmedians=True, widths=0.8)
        ax.set_xticklabels(labels, fontsize=size_font-5)
        violin['cmeans'].set_color('r')
        violin['cmeans'].set_linewidths(2)
        violin['cmedians'].set_color('g')
        for i in range(len(violin['bodies'])):
            violin['bodies'][i].set_color(set_color(i, 5))
            violin['bodies'][i].set_alpha(0.3)
            violin['bodies'][i].set_edgecolor(set_color(i, 5))

    elif down_dict is not None and up_down_dict is None:
        labels = ["up", "down", "ACE", "CCE"] * 4
        pos = [1, 2, 3, 4, 6, 7, 8, 9, 11, 12, 13, 14, 16, 17, 18, 19]
        ax.set_xticks(pos)
        box = ax.boxplot([up_dict[A], down_dict[A], dic_ACE[A], dic_CCE[A],
                          up_dict[C], down_dict[C], dic_ACE[C], dic_CCE[C],
                          up_dict[G], down_dict[G], dic_ACE[G], dic_CCE[G],
                          up_dict[T], down_dict[T], dic_ACE[T], dic_CCE[T]], positions=pos, labels=labels,
                         showmeans=True, notch=False, patch_artist=True)
        for i in range(len(box['boxes'])):
            box['boxes'][i].set_facecolor(set_color(i, 4))
        violin = ax.violinplot([up_dict[A], down_dict[A], dic_ACE[A], dic_CCE[A],
                                up_dict[C], down_dict[C], dic_ACE[C], dic_CCE[C],
                                up_dict[G], down_dict[G], dic_ACE[G], dic_CCE[G],
                                up_dict[T], down_dict[T], dic_ACE[T], dic_CCE[T]],
                               pos, showmeans=True, showmedians=True, widths=0.8)
        ax.set_xticklabels(labels, fontsize=size_font-5)
        violin['cmeans'].set_color('r')
        violin['cmeans'].set_linewidths(2)
        violin['cmedians'].set_color('g')
        for i in range(len(violin['bodies'])):
            violin['bodies'][i].set_color(set_color(i, 4))
            violin['bodies'][i].set_alpha(0.3)
            violin['bodies'][i].set_edgecolor(set_color(i, 4))

    else:
        labels = ["input", "ACE", "CCE"] * 4
        pos = [1, 2, 3, 5, 6, 7, 9, 10, 11, 13, 14, 15]
        ax.set_xticks(pos)
        box = ax.boxplot([up_dict[A], dic_ACE[A], dic_CCE[A],
                          up_dict[C], dic_ACE[C], dic_CCE[C],
                          up_dict[G], dic_ACE[G], dic_CCE[G],
                          up_dict[T], dic_ACE[T], dic_CCE[T]],
                         positions=pos, labels=labels, showmeans=True, notch=False, patch_artist=True)
        for i in range(len(box['boxes'])):
            box['boxes'][i].set_facecolor(set_color(i, 3))
        violin = ax.violinplot([up_dict[A], dic_ACE[A], dic_CCE[A],
                                up_dict[C], dic_ACE[C], dic_CCE[C],
                                up_dict[G], dic_ACE[G], dic_CCE[G],
                                up_dict[T], dic_ACE[T], dic_CCE[T]],
                               pos, showmeans=True, showmedians=True, widths=0.8)
        violin['cmeans'].set_color('r')
        violin['cmeans'].set_linewidths(2)
        violin['cmedians'].set_color('g')
        ax.set_xticklabels(labels)
        for i in range(len(violin['bodies'])):
            violin['bodies'][i].set_color(set_color(i, 3))
            violin['bodies'][i].set_alpha(0.3)
            violin['bodies'][i].set_edgecolor(set_color(i, 3))

    hb, = ax.plot([100, 100], 'b-')
    hr, = ax.plot([100, 100], 'r-')
    hy, = ax.plot([100, 100], 'y-')
    hg, = ax.plot([100, 100], 'g-')
    ax.plot([100, 100], 'w-', lw=2)

    plt.xlabel(u"Jeux d'exons", fontsize=size_font)
    if rich == "rich":
        ax.legend((hb, hr, hy, hg), ('A rich codons', 'C rich codons', 'G rich codons', 'T rich codons'))
        plt.title(u"Proportion de codons riches en A, C, G, T dans différents jeux d'exons", fontsize=size_font)
        plt.ylabel(u"Fréquences dans les exons", fontsize=size_font)
        plt.savefig(output + "proportion_nucleotide_rich_codon.png", bbox_inches='tight')
    elif rich == "poor":
        ax.legend((hb, hr, hy, hg), ('codon only B ', 'codon only D', 'codon only H', 'codon only V'))
        plt.title(u"Proportion de codons ne contenant que des nucleotides B, D, H, V"
                  u" dans différents jeux d'exons", fontsize=size_font)
        plt.ylabel(u"Fréquences dans les exons", fontsize=size_font)
        plt.savefig(output + "proportion_nucleotide_poor_codon.png", bbox_inches='tight')
    else:
        ax.legend((hb, hr, hy, hg), ('A ending codons', 'C ending codons', 'G ending codons', 'T ending codons'))
        plt.title(u"Proportion de codons finissant par A, C, G, T dans différents jeux d'exons", fontsize=size_font)
        plt.ylabel(u"Fréquences du dernier nucléotide dans les exons", fontsize=size_font)
        plt.savefig(output + "last_nucleotide_proportion.png", bbox_inches='tight')
    plt.clf()
    plt.cla()
    plt.close()


def boxplot_dinucleotide_rich_or_poor_codon(output, up_dict, down_dict, up_down_dict, rich_codon):
    """
    :param output: (string) the folder where the graphic will be created
    :param up_dict: (dictionary of list of float), a dictionary that gives for each up-regulated exons within a
    particular condition the proportion of codon rich/poor in dinucleotides.
    A codon rich in AG, for example, is a codon that contains at least 2 A or 2 G or AG.
    A codon poor in AG is a codon that contains no A and no G
    :param down_dict: (dictionary of list of float or None if not present), a dictionary that gives for each
    down-regulated exons within a particular condition the proportion of codon rich/poor in  dinucleotides.
    A codon rich in AG, for example, is a codon that contains at least 2 A or 2 G or AG.
    A codon poor in AG is a codon that contains no A and no G
    :param up_down_dict: (dictionary of list of float or None if not present), a dictionary that gives for each
    regulated (up+down) exons within a particular condition the proportion of codon rich/poor in dinucleotides.
    A codon rich in AG, for example, is a codon that contains at least 2 A or 2 G or AG.
    A codon poor in AG is a codon that contains no A and no G
    :param rich_codon: (string) "rich" if we want to create a graphics that gives the proportion of codons codons rich
    in dinuleotides in the up-down-up_and_down-regulated set of exons. "poor" if we want to create a graphics that gives
    the proportion of codons that contains no A, C, G, or T in the up-down-up_and_down-regulated set of exons. "last" if
    we want to create a graphics that gives the proportion of codons that ends with A, C, G, or T in the up-down-up_and_
    down-regulated set of exons
    """
    size_font = 20
    fig, ax = plt.subplots(figsize=(48. / 2.54, 27 / 2.54))
    # loading the variable we need to complete the graphics
    file_dir = os.path.dirname(os.path.realpath(__file__))
    sys.path.insert(0, file_dir + "/nucleotides_riche_codon_ACE_CCE/")
    if rich_codon == "rich":
        ACE_f = __import__("ACE_rich")
        CCE_f = __import__("CCE_rich")
    else:
        ACE_f = __import__("ACE_only")
        CCE_f = __import__("CCE_only")
    dic_ACE = ACE_f.dic_ACE
    dic_CCE = CCE_f.dic_CCE

    if down_dict is not None and up_down_dict is not None:
        labels = ["up", "down", "up\nand\ndown", "ACE", "CCE"] * 6
        pos = [1, 2, 3, 4, 5, 7, 8, 9, 10, 11, 13, 14, 15, 16, 17, 19, 20, 21, 22, 23, 25, 26, 27, 28, 29,
                       31, 32, 33, 34, 35]
        ax.set_xticks(pos)
        box = ax.boxplot([up_dict['Y'], down_dict['Y'], up_down_dict['Y'], dic_ACE['Y'], dic_CCE['Y'],
                          up_dict['R'], down_dict['R'], up_down_dict['R'], dic_ACE['R'], dic_CCE['R'],
                          up_dict['W'], down_dict['W'], up_down_dict['W'], dic_ACE['W'], dic_CCE['W'],
                          up_dict['S'], down_dict['S'], up_down_dict['S'], dic_ACE['S'], dic_CCE['S'],
                          up_dict['K'], down_dict['K'], up_down_dict['K'], dic_ACE['K'], dic_CCE['K'],
                          up_dict['M'], down_dict['M'], up_down_dict['M'], dic_ACE['M'], dic_CCE['M'],
                          ], positions=pos, labels=labels, showmeans=True, notch=False, patch_artist=True)
        for i in range(len(box['boxes'])):
            box['boxes'][i].set_facecolor(set_color(i, 5))
        violin = ax.violinplot([up_dict['Y'], down_dict['Y'], up_down_dict['Y'], dic_ACE['Y'], dic_CCE['Y'],
                                up_dict['R'], down_dict['R'], up_down_dict['R'], dic_ACE['R'], dic_CCE['R'],
                                up_dict['W'], down_dict['W'], up_down_dict['W'], dic_ACE['W'], dic_CCE['W'],
                                up_dict['S'], down_dict['S'], up_down_dict['S'], dic_ACE['S'], dic_CCE['S'],
                                up_dict['K'], down_dict['K'], up_down_dict['K'], dic_ACE['K'], dic_CCE['K'],
                                up_dict['M'], down_dict['M'], up_down_dict['M'], dic_ACE['M'], dic_CCE['M']],
                               pos, showmeans=True, showmedians=True, widths=0.8)
        violin['cmeans'].set_color('r')
        violin['cmeans'].set_linewidths(2)
        violin['cmedians'].set_color('g')
        ax.set_xticklabels(labels, fontsize=size_font-8)
        for i in range(len(violin['bodies'])):
            violin['bodies'][i].set_color(set_color(i, 5))
            violin['bodies'][i].set_alpha(0.3)
            violin['bodies'][i].set_edgecolor(set_color(i, 5))

    elif down_dict is not None and up_down_dict is None:
        labels = ["up", "down", "ACE", "CCE"] * 6
        pos = [1, 2, 3, 4, 6, 7, 8, 9, 11, 12, 13, 14, 16, 17, 18, 19, 21, 22, 23, 24, 26, 27, 28, 29]
        ax.set_xticks(pos)
        box = ax.boxplot([up_dict['Y'], down_dict['Y'], dic_ACE['Y'], dic_CCE['Y'],
                          up_dict['R'], down_dict['R'], dic_ACE['R'], dic_CCE['R'],
                          up_dict['W'], down_dict['W'], dic_ACE['W'], dic_CCE['W'],
                          up_dict['S'], down_dict['S'], dic_ACE['S'], dic_CCE['S'],
                          up_dict['K'], down_dict['K'], dic_ACE['K'], dic_CCE['K'],
                          up_dict['M'], down_dict['M'], dic_ACE['M'], dic_CCE['M'],
                          ], positions=pos, labels=labels, showmeans=True, notch=False, patch_artist=True)
        for i in range(len(box['boxes'])):
            box['boxes'][i].set_facecolor(set_color(i, 4))
        violin = ax.violinplot([up_dict['Y'], down_dict['Y'], dic_ACE['Y'], dic_CCE['Y'],
                                up_dict['R'], down_dict['R'], dic_ACE['R'], dic_CCE['R'],
                                up_dict['W'], down_dict['W'], dic_ACE['W'], dic_CCE['W'],
                                up_dict['S'], down_dict['S'], dic_ACE['S'], dic_CCE['S'],
                                up_dict['K'], down_dict['K'], dic_ACE['K'], dic_CCE['K'],
                                up_dict['M'], down_dict['M'], dic_ACE['M'], dic_CCE['M']],
                               pos, showmeans=True, showmedians=True, widths=0.8)
        violin['cmeans'].set_color('r')
        violin['cmeans'].set_linewidths(2)
        violin['cmedians'].set_color('g')
        ax.set_xticklabels(labels, fontsize=size_font-8)
        for i in range(len(violin['bodies'])):
            violin['bodies'][i].set_color(set_color(i, 4))
            violin['bodies'][i].set_alpha(0.3)
            violin['bodies'][i].set_edgecolor(set_color(i, 4))

    else:
        labels = ["input", "ACE", "CCE"] * 6
        pos = [1, 2, 3, 5, 6, 7, 9, 10, 11, 13, 14, 15, 17, 18, 19, 21, 22, 23]
        ax.set_xticks(pos)
        box = ax.boxplot([up_dict['Y'], dic_ACE['Y'], dic_CCE['Y'],
                          up_dict['R'], dic_ACE['R'], dic_CCE['R'],
                          up_dict['W'], dic_ACE['W'], dic_CCE['W'],
                          up_dict['S'], dic_ACE['S'], dic_CCE['S'],
                          up_dict['K'], dic_ACE['K'], dic_CCE['K'],
                          up_dict['M'], dic_ACE['M'], dic_CCE['M']],
                         positions=pos, labels=labels, showmeans=True, notch=False, patch_artist=True)
        for i in range(len(box['boxes'])):
            box['boxes'][i].set_facecolor(set_color(i, 3))
        violin = ax.violinplot([up_dict['Y'], dic_ACE['Y'], dic_CCE['Y'],
                                up_dict['R'], dic_ACE['R'], dic_CCE['R'],
                                up_dict['W'], dic_ACE['W'], dic_CCE['W'],
                                up_dict['S'], dic_ACE['S'], dic_CCE['S'],
                                up_dict['K'], dic_ACE['K'], dic_CCE['K'],
                                up_dict['M'], dic_ACE['M'], dic_CCE['M']],
                               pos, showmeans=True, showmedians=True, widths=0.8)
        violin['cmeans'].set_color('r')
        violin['cmeans'].set_linewidths(2)
        violin['cmedians'].set_color('g')
        ax.set_xticklabels(labels, fontsize=size_font-8)
        for i in range(len(violin['bodies'])):
            violin['bodies'][i].set_color(set_color(i, 3))
            violin['bodies'][i].set_alpha(0.3)
            violin['bodies'][i].set_edgecolor(set_color(i, 3))
    hb, = ax.plot([100, 100], 'b-')
    hr, = ax.plot([100, 100], 'r-')
    hy, = ax.plot([100, 100], 'y-')
    hg, = ax.plot([100, 100], 'g-')
    hk, = ax.plot([100, 100], 'k-')
    hc, = ax.plot([100, 100], 'c-')
    ax.plot([100, 100], 'w-', lw=2)

    plt.xlabel(u"Jeux d'exons", fontsize=size_font)
    plt.ylabel(u"Fréquences dans les exons", fontsize=size_font)
    if rich_codon == "rich":
        ax.legend((hb, hr, hy, hg, hk, hc), ('Y rich codons', 'R rich codons', 'W rich codons', 'S rich codons',
                                             'K rich codon', 'M rich codon'), loc=1)
        plt.title(u"Proportion en codons riches en nucleotides Y, R, W, S, K, M dans différents jeux d'exons",
                  fontsize=size_font)
        plt.savefig(output + "proportion_dinucleotide_rich_codon.png", bbox_inches='tight')
    else:
        ax.legend((hb, hr, hy, hg, hk, hc), ('codons only Y', 'codons only R', 'codons only W', 'codon only S',
                                             'codon only K', 'codons only M'), loc=1)
        plt.title(u"Proportion en codons ne contenant que des nucléotides  Y, R, W, S, K, M dans différents jeux "
                  u"d'exons",
                  fontsize=size_font)
        plt.savefig(output + "proportion_dinucleotide_only_codon.png", bbox_inches='tight')

    plt.clf()
    plt.cla()
    plt.close()
