#!/usr/bin/env python
# coding: utf-8

# for each hla, whether the TCRs from all positive pairs appears in
# at least two subjects with known exist or not status for that hla


import os
import sys

import math

import random

import argparse

import numpy as np
import pandas as pd

from collections import defaultdict

sys.path.append(os.getcwd())


parser = argparse.ArgumentParser(description='sample negative TCRs for each HLA allele')
parser.add_argument('--hla_class', default="HLA_I", help='the HLA class to use')
parser.add_argument('--hla_i', type=int, default=1, help='hla allele index, staring from 1')


def sample_tcrs(hla_class="HLA_I", hla_i=1):

    hla_class = args.hla_class
    hla_i = args.hla_i
    hla_i = int(hla_i)

    print(hla_class)
    print(hla_i)

    n_fold = 5

    data_dir = "../../../data/intermediate_data/"

    if hla_class == "HLA_I":
        asso_file = data_dir + "t12_HLA_I_associated_TCR_v_alleles.csv"
    elif hla_class == "HLA_II":
        asso_file = data_dir + "t12_HLA_II_no_haplotype_associated_TCR_v_alleles.csv"
    else:
        print("hla_class input wrong. should be either HLA_I or HLA_II")
        return

    df_asso = pd.read_csv(asso_file, sep=",", header=0)

    # create a dict to hold the set of TCRs that each HLA is associated with

    pos_set_dict = defaultdict(set)

    for a, b in zip(df_asso.tcr.tolist(), df_asso.hla_allele.tolist()):
        pos_set_dict[b].add(a)

    # identify HLA alleles under the order in orginal file

    HLA_v2_features = pd.read_csv(data_dir + "HLA_v2_features_reformat.csv", header=0)

    if hla_class == "HLA_I":
        HLA_I_short = set(["HLA-A", "HLA-B", "HLA-C"])
        HLA_list = [x for x in HLA_v2_features.hla.tolist() if x[:5] in HLA_I_short]
    else:
        haplotypes = ['HLA-DRDQ*15:01_01:02_06:02', 'HLA-DRDQ*03:01_05:01_02:01',
                      'HLA-DRDQ*13:01_01:03_06:03', 'HLA-DRDQ*10:01_01:05_05:01',
                      'HLA-DRDQ*09:01_03:02_03:03']
        HLA_list = [x for x in HLA_v2_features.hla.tolist() if ((x[:5] == "HLA-D") and (x not in haplotypes))]

    cur_hla = HLA_list[hla_i-1]

    pos_tcr_list = list(pos_set_dict[cur_hla])

    # get the frequency of positive TCRs among subjects with known status for each HLA

    df_tcr_freq = pd.read_csv(data_dir + hla_class + "_all_match_prepare/tcr_freq_" + hla_class + "_"+str(hla_i)+".csv",
                              sep=",", header=0)

    # build a dictionary with full TCR name as key and row index as value
    row_id_dict = defaultdict(int)

    for row_i, (v, aa) in enumerate(zip(df_tcr_freq.v_allele, df_tcr_freq.amino_acids)):
        tcr_full = v + "," + aa
        row_id_dict[tcr_full] = row_i

    greaterthan2flag = sum([x in row_id_dict for x in pos_tcr_list])/len(pos_tcr_list)

    df_flag = pd.DataFrame([greaterthan2flag], columns=["flag"])
    df_flag.to_csv(data_dir + hla_class + "_all_match_prepare_2_nolessthan2_flags/" + "neg_pairs_"+str(hla_i)+"_flag.csv",
                   index=False)

if __name__ == "__main__":
    args = parser.parse_args()
    sample_tcrs(args)
