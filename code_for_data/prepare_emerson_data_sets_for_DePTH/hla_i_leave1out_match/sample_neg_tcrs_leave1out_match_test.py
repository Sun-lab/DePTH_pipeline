#!/usr/bin/env python
# coding: utf-8

# overall goal ------------------------
# the goal is to put only one HLA allele in testing data, while make sure that
# the tcrs in training + validation do not have overlap with those in test
# still want to sample according to frequency matching

# put one target HLA allele in test data,
# for training and vaidation positive pairs, only keep those involving TCRs not
# appearing in positive test (why the train + valid pos pairs are decided first,
# is because the number of positive pairs is very limited)
# construct test negative pairs, such that the tcrs do not overlap with whose
# in training + validation positive
# for training and validatoin negative paris, construct in a way that do not
# involve any tcrs in either test positive or testing negative

# for each HLA
# separate the tcrs associated with given HLA using percentiles of frequency
# sample n_fold * matching tcrs from the negative
# with exceptions for extremely low frequency, 0, 1, 2, 3, up to 10, and then
# folowing percentiles at 10%, 20%, 30%, ...
# if there are not enough TCRs under a certain interval, just reuse the
# negative tcrs

# goal for this file -----------------------
# output test positive/negative pairs
# output train_valid positive pairs
# output a correspondance data frame of hla name and hla_i (1-indexed)

import os
import sys

import math

import random

import argparse

import numpy as np
import pandas as pd

from random import shuffle
from collections import defaultdict

sys.path.append(os.getcwd())

parser = argparse.ArgumentParser(description='sample negative TCRs for each HLA-I allele')
parser.add_argument('--hla_class', default='HLA_I', help='the HLA class to use')
parser.add_argument('--test_hla', default='HLA-B_08_01',
                    help='the HLA allele to leave out in testing data')
parser.add_argument('--label', default='leave1out_match',
                    help='part of the label of the model to tell different settings apart')
parser.add_argument('--n_fold', default=5, type=int,
                     help='number of negative pairs divided by that of positive pairs')


def prepare_for_sample(hla_class='HLA_I', test_hla='HLA-B_08_01',
                       label='leave1out_match', n_fold=5):
    hla_class = args.hla_class
    test_hla = args.test_hla
    label = args.label
    n_fold = int(args.n_fold)

    test_hla_strs = test_hla.split("_")
    test_hla = test_hla_strs[0] + "*" + test_hla_strs[1] + ":" + test_hla_strs[2]

    name = hla_class + "_" + label

    data_dir = "../../../data/intermediate_data/"

    # load the association matrix with specified hla class
    if hla_class == "HLA_I":
        asso_file = data_dir + "t12_HLA_I_associated_TCR_v_alleles.csv"
    elif hla_class == "HLA_II":
        asso_file = data_dir + "t12_HLA_II_no_haplotype_associated_TCR_v_alleles.csv"
    else:
        print("hla_class input wrong. should be either HLA_I or HLA_II")
        return

    df_asso = pd.read_csv(asso_file, sep=",", header=0)

    # positive pair list:
    positive_ori = [(tcr, hla) for tcr, hla in zip(df_asso.tcr.to_list(), df_asso.hla_allele.to_list())]
    # positives involving the HLA allele to leave out in test data
    positive_test_ori = [pair for pair in positive_ori if pair[1] == test_hla]
    # tcrs that we want to exclude from training and validation positive pairs
    tcr_to_exclude = set([pair[0] for pair in positive_test_ori])
    positive_left_ori = [pair for pair in positive_ori if
                         ((pair[1] != test_hla) and (pair[0] not in tcr_to_exclude))]

    # identify HLA-I alleles under the order in orginal file

    HLA_v2_features = pd.read_csv(data_dir + "HLA_v2_features_reformat.csv",
                                  header=0)
    HLA_I_short = set(["HLA-A", "HLA-B", "HLA-C"])
    if hla_class == "HLA_I":
        HLA_list = [x for x in HLA_v2_features.hla.tolist() if x[:5] in HLA_I_short]
    elif hla_class == "HLA_II":
        HLA_list = [x for x in HLA_v2_features.hla.tolist() if x[:5] not in HLA_I_short]
    else:
        print("HLA class not coded for yet")
        return
    # find the index of the test hla
    # to keep it consistent with other code files, make the hla_i 1-indexed
    test_hla_i = [i for i, hla in enumerate(HLA_list) if hla == test_hla][0] + 1

    df_hla_index = pd.DataFrame(zip(HLA_list, list(range(1, len(HLA_list) + 1))),
                                columns=['hla', 'index1'])

    # prepare folder for writing files out to
    output_dir = data_dir + hla_class + "_leave1out_match_prepare"

    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)

    df_hla_index.to_csv(output_dir + "/" + hla_class + "_index1.csv", index=False)

    df_pos_left_ori = pd.DataFrame(positive_left_ori, columns=['tcr', 'hla_allele'])
    df_pos_left_ori.to_csv(output_dir + "/" + "pos_left_ori_" + str(test_hla_i) + ".csv", index=False)
    # convert the set of tcrs from positive test to list
    # to prepare for frequency matching
    pos_tcr_list = list(tcr_to_exclude)

    # get the frequency of positive TCRs among subjects with known status
    # for the current test hla
    df_tcr_freq = pd.read_csv(data_dir + hla_class + "_all_match_prepare/tcr_freq_"
                              + hla_class + "_" + str(test_hla_i) + ".csv",
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

    # all the tcrs involved in positive pairs need to be excluded from the
    # testing negative pairs
    pos_all_row_ids = set([row_id_dict[x] for x in df_asso.tcr.to_list()])

    # get the list of frequency for the tcrs involved in positive test pairs

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

    for bin_id in range(len(all_grids) - 1):
        left_freq = all_grids[bin_id]
        right_freq = all_grids[bin_id + 1]
        bin_list = [i for i, frequency in enumerate(df_tcr_freq.freq) if
                    ((frequency > left_freq) and (frequency <= right_freq)
                     and (i not in pos_all_row_ids))]
        neg_bin_dict[str(bin_id)] = bin_list
        pos_bin_len_dict[str(bin_id)] = sum([1 for frequency in pos_freq if
                                             (frequency > left_freq and frequency <= right_freq)])

    # load random seeds, one for each HLA
    df_seeds = pd.read_csv(data_dir + hla_class + "_all_match_random_seeds.csv",
                           header=None)
    cur_seed = df_seeds[0][test_hla_i - 1]

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
            n_repeat = math.ceil(bin_num_neg / len(bin_tcr_list))
            bin_tcr_list_to_choose = [x for _ in range(n_repeat) for x in bin_tcr_list]
        else:
            bin_tcr_list_to_choose = bin_tcr_list
        neg_row_ids_list += bin_tcr_list_to_choose[:bin_num_neg]

    neg_v_allele_list = [df_tcr_freq.v_allele[i] for i in neg_row_ids_list]
    neg_amino_acids_list = [df_tcr_freq.amino_acids[i] for i in neg_row_ids_list]

    neg_tcr_full_list = \
        [v + "," + aa for v, aa in zip(neg_v_allele_list, neg_amino_acids_list)]

    cur_hla_list = [test_hla for _ in range(len(neg_row_ids_list))]

    shuffle(neg_tcr_full_list)
    shuffle(cur_hla_list)
    shuffle(positive_test_ori)

    df_test_neg = pd.DataFrame(list(zip(neg_tcr_full_list, cur_hla_list)),
                               columns=["tcr", "hla_allele"])
    df_test_pos = pd.DataFrame(positive_test_ori, columns=['tcr', 'hla_allele'])

    model_data_dir = data_dir + name + "_raw_v"
    if not os.path.isdir(model_data_dir):
        os.mkdir(model_data_dir)

    df_test_pos.to_csv(model_data_dir + "/" + name + "_test_pos_" + \
                       str(test_hla_i) + ".csv", index=False)
    df_test_neg.to_csv(model_data_dir + "/" + name + "_test_neg_" + \
                       str(test_hla_i) + ".csv", index=False)

    df_flag = pd.DataFrame([insufficient_tcr_flag], columns=["flag"])
    df_flag.to_csv(output_dir + "/" +
                   "neg_pairs_" + str(test_hla_i) + "_flag.csv",
                   index=False)


if __name__ == "__main__":
    args = parser.parse_args()
    prepare_for_sample(args)
