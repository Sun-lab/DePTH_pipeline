#!/usr/bin/env python
# coding: utf-8

# separate the tcrs associated with given HLA using percentiles of frequency
# sample n_fold * matching tcrs from the negative
# with exceptions for extremely low frequency, 0, 1, 2, 3, up to 10, and then
# folowing percentiles at 10%, 20%, 30%, ...
# if there are not enough TCRs under a certain interval, just reuse the
# negative tcrs


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

    # build a dictionary with row index as key and freq as value
    freq_dict = defaultdict(int)

    for row_i, frequency in enumerate(df_tcr_freq.freq):
        freq_dict[row_i] = frequency

    # it has been verified in a previous step that for each HLA allele,
    # the tcrs from the associated pairs each appears in at least two subjects
    # with known existing or not status for this HLA

    # get the set of row index for the tcrs involved in positive pairs
    # get the list of frequency for the tcrs involved in positive pairs

    pos_row_ids = set([row_id_dict[x] for x in pos_tcr_list])
    pos_freq = [freq_dict[row_id_dict[x]] for x in pos_tcr_list]

    np_pos_freq = np.array(pos_freq)

    percentile_list = []
    percentile_list += [0]
    for percent_id in range(1, 10):
        percentile_list += [np.percentile(np_pos_freq, 10 * percent_id, interpolation='linear')]

    percentile_list += [max(np_pos_freq)]

    # get the first batch of splits
    # there is no tcr with 0 frequency, and most likely no with 1 frequency
    # but still keep it here, since it should not have effect on the result
    first_splits = list(range(11))
    # find the first number from percentile_list that the max of first splits
    # is smaller than
    ind_grt = sum([x <= max(first_splits) for x in percentile_list])

    all_grids = first_splits + percentile_list[ind_grt:]

    # get a dictionary to hold the number of positive tcrs
    # for each bin

    # get a dictionary to hold a list of tcrs for each bin
    # exclude those already appeared in positive

    pos_bin_len_dict = defaultdict(int)
    neg_bin_dict = defaultdict(list)

    for bin_id in range(len(all_grids)-1):
        left_freq = all_grids[bin_id]
        right_freq = all_grids[bin_id+1]
        bin_list = [i for i, frequency in enumerate(df_tcr_freq.freq) if
                    ((frequency > left_freq) and (frequency <= right_freq)
                     and (i not in pos_row_ids))]
        neg_bin_dict[str(bin_id)] = bin_list
        pos_bin_len_dict[str(bin_id)] = sum([1 for frequency in pos_freq if
                                            ((frequency > left_freq) and (frequency <= right_freq))])

    # load random seed file
    df_seeds = pd.read_csv(data_dir + hla_class + "_all_match_random_seeds.csv", header=None)
    cur_seed = df_seeds[0][hla_i-1]

    random.seed(int(cur_seed))

    # shuffle the row indexes in each bin
    # take the first certain number of them
    # it is possible that some bin may not have enough unique TCRs to offer
    # in this case, currently we reuse some TCRs

    neg_row_ids_list = []

    # use a flag to keep track of whether there is the situation of not
    # sufficiently many negative TCRs in the matched freq interval
    insufficient_tcr_flag = 0

    for str_bin_id in pos_bin_len_dict:
        bin_num_pos = pos_bin_len_dict[str_bin_id]
        bin_num_neg = n_fold * bin_num_pos
        bin_tcr_list = neg_bin_dict[str_bin_id]
        random.shuffle(bin_tcr_list)
        if bin_num_neg > len(bin_tcr_list):
            insufficient_tcr_flag = 1
            n_repeat = math.ceil(bin_num_neg/len(bin_tcr_list))
            bin_tcr_list_to_choose = [x for _ in range(n_repeat) for x in bin_tcr_list]
        else:
            bin_tcr_list_to_choose = bin_tcr_list
        neg_row_ids_list += bin_tcr_list_to_choose[:bin_num_neg]

    neg_v_allele_list = [df_tcr_freq.v_allele[i] for i in neg_row_ids_list]
    neg_amino_acids_list = [df_tcr_freq.amino_acids[i] for i in neg_row_ids_list]

    neg_tcr_full_list = \
        [v+","+aa for v, aa in zip(neg_v_allele_list, neg_amino_acids_list)]

    cur_hla_list = [cur_hla for _ in range(len(neg_row_ids_list))]

    df_output = pd.DataFrame(list(zip(neg_tcr_full_list, cur_hla_list)),
                             columns=["tcr", "hla"])

    output_fname = "neg_pairs_"+str(hla_i)+".csv"
    df_output.to_csv(data_dir + hla_class + "_all_match_prepare_2/"+output_fname,
                     index=False)

    df_flag = pd.DataFrame([insufficient_tcr_flag], columns=["flag"])
    df_flag.to_csv(data_dir + hla_class + "_all_match_prepare_2_flags/" + "neg_pairs_"+str(hla_i)+"_flag.csv",
                   index=False)

if __name__ == "__main__":
    args = parser.parse_args()
    sample_tcrs(args)
