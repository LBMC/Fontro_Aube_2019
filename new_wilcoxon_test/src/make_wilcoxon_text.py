#!/usr/bin/python3

"""
Description:
    The goal of this script is to make the wilcoxon test on the new figures.
"""

import pandas as pd
import rpy2.robjects as robj
import rpy2.robjects.vectors as v
import subprocess
import os
import numpy as np


def get_values_from_files(enrichment_file, sheet, row_name):
    """
    Get the interest value and the control value in ``enrichment_file`` stored in the sheet ``sheet`` under the \
    row named ``row_named`` and the columns "frequencies_of_the_interest_set" and \
    "average_frequencies_of_the_10000_sets" respectively

    :param enrichment_file: (string) a tRNA enrichment file
    :param sheet: (string) the name of the wanted sheet
    :param row_name: (string) the name of the interest row
    :return: (2 int) the test and the control
    """
    xl = pd.ExcelFile(enrichment_file)
    # opening the sheet of interest
    df = xl.parse(sheet)
    test_val = df.loc[df[sheet.replace("nt_pos", "ntp")] == row_name, ["frequencies_of_the_interest_set",
                                                                       "average_frequencies_of_the_10000_sets"]].head(1)
    return float(test_val.values[0][0]), float(test_val.values[0][1])


def mann_withney_test_r(list_values1, list_values2):
    wicox = robj.r("""

    function(x, y){
        test = wilcox.test(x,y, alternative='two.sided', correct=F)
        return(test$p.value)
    }

                   """)
    pval = float(wicox(v.FloatVector(list_values1), v.FloatVector(list_values2))[0])
    return pval


def get_list_val(list_files, sheet, row_name):
    """
    Get the list of target and control values of interest.

    :param list_files: (list of string) list of tRNA enrichment file
    :param sheet: (string) the sheet of interest
    :param row_name: (string) the target row wanted
    :return: (2 list of values) list of test values and list of control values
    """
    test_list = []
    ctrl_list = []
    for enrichment_file in list_files:
        test, ctrl = get_values_from_files(enrichment_file, sheet, row_name)
        test_list.append(test)
        ctrl_list.append(ctrl)
    return test_list, ctrl_list


def get_file_list(splicing_factor, regulation):
    """
    Get the enrichment files corresponding to the splicing factor ``splicing_factor``.

    :param splicing_factor: (string) a splicing factor
    :param regulation: (string) up or down
    :return: (list of string) the list of files wanted
    """
    target_path = os.path.realpath(os.path.dirname(__file__)).replace("src", "data/")
    find_path = target_path + splicing_factor
    res_file = subprocess.check_output("find %s -name enrichment_report.xlsx -type f" % find_path, shell=True,
                                       stderr=subprocess.STDOUT).decode("ascii").strip("\n").split("\n")
    new_res = []
    list_name = []
    for my_file in res_file:
        if regulation in my_file:
            new_res.append(my_file)
            list_name.append("_".join(my_file.split("/")[-4].split("_")[2:]))
    return new_res, list_name


def main():
    analysis = {"SRSF1": {"sheet": ["nt_info", "nt_pos_info", "nt_pos_info", "amino_acid",
                                    "amino_acid", "amino_acid", "amino_acid"],
                          "row": ["G", "G3", "G1n2", "G", "A", "P", "K"],
                          "fig": ["Fig_1A", "Fig_2B_G3", "Fig_2B_G1-2", "Fig_3B_G",
                                  "Fig_3B_A", "Fig_3B_P", "Fig_3B_K"],
                          "regulation": ["down", "up"]},
                "hnRNPL": {"sheet": ["feature", "feature"],
                           "row": ["Hydroxylic", "Negatively-charged"],
                           "fig": ["Fig_5F_Hydroxil_hnRNPL", "Fig_5F_Negatively_charged_hnRNPL"],
                           "regulation": ["up"]},
                "PTBP1": {"sheet": ["feature", "feature"],
                          "row": ["Hydroxylic", "Negatively-charged"],
                          "fig": ["Fig_5F_Hydroxil_PTBP1", "Fig_5F_Negatively_charged_PTBP1"],
                          "regulation": ["up"]}
                }
    res_dic = {}
    output = os.path.realpath(os.path.dirname(__file__)).replace("src", "result/fig_stat.txt")
    outfile = open(output, "w")
    for splicing_factor in analysis.keys():
        for regulation in analysis[splicing_factor]["regulation"]:
            for i in range(len(analysis[splicing_factor]["sheet"])):
                print("%s %s %s" % (splicing_factor, regulation, analysis[splicing_factor]["sheet"][i]))
                list_files, name_files = get_file_list(splicing_factor, regulation)
                target_val, control_val = get_list_val(list_files, analysis[splicing_factor]["sheet"][i],
                                                       analysis[splicing_factor]["row"][i])
                relative = (np.array(target_val) - np.array(control_val)) / np.array(control_val) * 100
                df = pd.DataFrame({"target": target_val, "ctrl": control_val, "relative": relative}, index=name_files)
                df = df.transpose().loc[["target", "ctrl", "relative"]]
                pval = mann_withney_test_r(target_val, control_val)
                key = "%s_%s" % (analysis[splicing_factor]["fig"][i], regulation)
                res_dic[key] = "-----------%s %s exons---------\n" % (analysis[splicing_factor]["fig"][i], regulation)
                res_dic[key] += df.to_string() + "\n"
                res_dic[key] += "pval = %s\n\n" % pval
    for fig in sorted(res_dic.keys()):
        outfile.write(res_dic[fig])
    outfile.close()


if __name__ == "__main__":
    main()
