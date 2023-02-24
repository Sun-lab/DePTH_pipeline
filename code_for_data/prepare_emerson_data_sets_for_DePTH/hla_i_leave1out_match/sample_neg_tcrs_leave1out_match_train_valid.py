#!/usr/bin/env python
# coding: utf-8


# overall goal ------------------------
# the goal is to put only one HLA allele in test data, while make sure that
# the tcrs in training + validation do not have overlap with those in test
# still want to sample according to frequency matching

# put one target HLA allele in test data,
# for training and vaidation positive pairs, only keep those involving TCRs not
# appearing in positive test (why the train + valid pos pairs are decided first,
# is because the number of positive pairs is very limited)
# construct test negative pairs, such that the tcrs do not overlap with whose
# in training + validation positive
# for training and validatoin negative paris, construct in a way that do not
# involve any tcrs in either test positive or test negative

# for each HLA
# separate the tcrs associated with given HLA using percentiles of frequency
# sample n_fold * matching tcrs from the negative
# with exceptions for extremely low frequency, 0, 1, 2, 3, up to 10, and then
# folowing percentiles at 10%, 20%, 30%, ...
# if there are not enough TCRs under a certain interval, just reuse the
# negative tcrs

# goal for this file -----------------------
# output train/valid positive/negative pairs


import os
import sys
import re
import pickle

import math

import random

import argparse

import numpy as np
import pandas as pd

from collections import defaultdict
from collections import Counter

sys.path.append(os.getcwd())


parser = argparse.ArgumentParser(description='sample negative TCRs for each HLA-I allele in train/valid')
parser.add_argument('--hla_class', default='HLA_I', help='the HLA class to use')
parser.add_argument('--test_hla_i', type=int, default=1,
                    help='the index (starting from 1) of the HLA allele to leave out in testing data')
parser.add_argument('--hla_i', type=int, default=2,
                    help='the index (starting from 1, cannot be the same as test_hla_i) of the HLA allele to sample tcr for')
parser.add_argument('--label', default='leave1out_match',
                    help='part of the label of the model to tell different settings apart')
parser.add_argument('--n_fold', default=5, type=int,
                     help='number of negative pairs divided by that of positive pairs')



def sample_tcrs(hla_i):

    hla_class = args.hla_class
    test_hla_i = int(args.test_hla_i)
    hla_i = int(args.hla_i)
    label = args.label
    n_fold = int(args.n_fold)

    if test_hla_i == hla_i:
        print("the hla allele to sample tcrs for cannot be the one in test")
        return

    name = hla_class + "_" + label

    data_dir = "../../../data/intermediate_data/"

    # load the corresponding file of positive training/validation pairs

    output_dir = data_dir + name + "_prepare"

    df_pos_left_ori = pd.read_csv(output_dir + "/" + "pos_left_ori_" + \
                                  str(test_hla_i) + ".csv", header=0)

    # create a dict to hold the set of TCRs that each HLA is associated with

    pos_set_dict = defaultdict(set)

    for a, b in zip(df_pos_left_ori.tcr.tolist(),
                    df_pos_left_ori.hla_allele.tolist()):
        pos_set_dict[b].add(a)

    # load the corresponding pairs for test
    model_data_dir = data_dir + name

    df_test_pos = pd.read_csv(model_data_dir + "/" + name + "_test_pos_" + \
                              str(test_hla_i) + ".csv", header=0)
    df_test_neg = pd.read_csv(model_data_dir + "/" + name + "_test_neg_" + \
                              str(test_hla_i) + ".csv", header=0)
    # this is the set of tcrs that cannot appear among the tcrs for train/valid
    tcrs_to_exclude = list(set(df_test_pos.tcr.tolist() + df_test_neg.tcr.tolist()))

    # load the file for HLA index, to identify the current HLA and associated TCRs

    df_hla_index = pd.read_csv(output_dir + "/" + hla_class + "_index1.csv",
                               header=0)

    cur_hla = df_hla_index.hla.tolist()[hla_i-1]

    pos_tcr_list = list(pos_set_dict[cur_hla])

    # get the frequency of positive TCRs among subjects with known status for each HLA

    df_tcr_freq = pd.read_csv(data_dir + hla_class + "_all_match_prepare/tcr_freq_"
                              + hla_class + "_" + str(hla_i) + ".csv",
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

    # get the set of row index for the tcrs involved in test pairs
    # get the set of row index for the tcrs involved in positive pairs
    # get the list of frequency for the tcrs involved in positive pairs

    tcrs2exclude_row_ids = set([row_id_dict[x] for x in tcrs_to_exclude])

    pos_row_ids = set([row_id_dict[x] for x in pos_tcr_list])
    pos_freq = [freq_dict[row_id_dict[x]] for x in pos_tcr_list]

    np_pos_freq = np.array(pos_freq)

    percentile_list = []
    percentile_list += [0]
    for percent_id in range(1,10):
        percentile_list += \
        [np.percentile(np_pos_freq, 10 * percent_id ,interpolation='linear')]

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
    # exclude those already appearing in positive
    # and those already appearing in test pos + neg

    pos_bin_len_dict = defaultdict(int)
    neg_bin_dict = defaultdict(list)

    for bin_id in range(len(all_grids)-1):
        left_freq = all_grids[bin_id]
        right_freq = all_grids[bin_id+1]
        bin_list = [i for i, frequency in enumerate(df_tcr_freq.freq) if
                    ((frequency > left_freq) and (frequency <= right_freq)
                     and (i not in pos_row_ids)
                     and (i not in tcrs2exclude_row_ids))]
        neg_bin_dict[str(bin_id)] = bin_list
        pos_bin_len_dict[str(bin_id)] = sum([1 for frequency in pos_freq if
                                            (frequency > left_freq and frequency <= right_freq)])


    # load random seeds, one for each HLA
    df_seeds = pd.read_csv(data_dir + hla_class + "_all_match_random_seeds.csv",
                           header=None)
    cur_seed = df_seeds[0][hla_i - 1]

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
            bin_tcr_list_to_choose = \
            [x for _ in range(n_repeat) for x in bin_tcr_list]
        else:
            bin_tcr_list_to_choose = bin_tcr_list
        neg_row_ids_list += bin_tcr_list_to_choose[:bin_num_neg]



    neg_v_allele_list = [df_tcr_freq.v_allele[i] for i in neg_row_ids_list]
    neg_amino_acids_list = [df_tcr_freq.amino_acids[i] for i in neg_row_ids_list]


    neg_tcr_full_list = \
        [v+","+aa for v, aa in zip(neg_v_allele_list, neg_amino_acids_list)]

    cur_hla_list = [cur_hla for _ in range(len(neg_row_ids_list))]



    df_output = pd.DataFrame(list(zip(neg_tcr_full_list, cur_hla_list)),
                             columns = ["tcr", "hla"])

    data_subfolder = data_dir + name + "_prepare_2/neg_pairs_test_hla_" + \
                       str(test_hla_i)

    if not os.path.isdir(data_subfolder):
        os.mkdir(data_subfolder)

    output_fname = "neg_pairs_" + str(hla_i) + ".csv"
    df_output.to_csv(data_subfolder + "/" + output_fname, index = False)

    flag_subfolder = data_dir + name + "_prepare_2_flags/neg_pairs_test_hla_" + \
                       str(test_hla_i)

    if not os.path.isdir(flag_subfolder):
        os.mkdir(flag_subfolder)

    df_flag = pd.DataFrame([insufficient_tcr_flag], columns = ["flag"])
    df_flag.to_csv(flag_subfolder + "/neg_pairs_"+str(hla_i)+"_flag.csv",
                   index = False)



if __name__ == "__main__":
    args = parser.parse_args()
    sample_tcrs(args)
