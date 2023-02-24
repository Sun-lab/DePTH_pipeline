#!/usr/bin/env python3

import os
import sys
import pandas as pd

import random
from random import shuffle

import argparse

# the following line is to deal with the case that self-defined module
# in the submission directory are not found when submitted through sbatch, as in
# https://stackoverflow.com/questions/46718499/
# how-does-one-make-sure-that-the-python-submission-script-in-slurm-is-in-the-loca/
# 46724189?noredirect=1#comment80425963_46724189
sys.path.append(os.getcwd())


parser = argparse.ArgumentParser(description='Generate (TCR, HLA) pairs for building the model')
parser.add_argument('--hla_class', default='HLA_I', help='the class of HLA')
parser.add_argument('--test_hla_i', type=int, default=1,
                    help='the index (starting from 1) of the HLA allele to leave out in testing data')
parser.add_argument('--limit_hla_i', type=int, default=85,
                    help='the max index (starting from 1) of the HLA alleles')
parser.add_argument('--label', default='leave1out_match',
                    help='part of the label of the model to tell different settings apart')
parser.add_argument('--prop_train', type=float, default=0.75,
                    help='proportion of train+valid pairs to be allocated to training')
parser.add_argument('--n_fold', type=int, default=5,
                    help='number of negative pairs divided by number of positive pairs')
parser.add_argument('--rseed', type=int, default=1627, help='random seed for random')


# generate pos/neg pairs for train/valid

def get_data(hla_class="HLA_I", test_hla_i=1, limit_hla_i=85,
             label="leave1out_match", prop_train=0.75, n_fold=5, rseed=1627):

    hla_class = args.hla_class
    test_hla_i = int(args.test_hla_i)
    limit_hla_i = int(args.limit_hla_i)
    label = args.label
    prop_train = float(args.prop_train)
    n_fold = int(args.n_fold)
    rseed = int(args.rseed)

    data_dir = "../../../data/intermediate_data/"
    # join hla class and part of the label to form a label for telling data resources apart
    name = hla_class + "_" + label

    # set random seed
    random.seed(rseed)

    # load the saved positive pairs for train + valid
    pos_file = \
        data_dir + name + "_prepare/pos_left_ori_" + str(test_hla_i) + ".csv"
    df_pos = pd.read_csv(pos_file, sep=",", header=0)

    positive_train_valid = [(tcr, hla) for tcr, hla in zip(df_pos.tcr.tolist(), df_pos.hla_allele.tolist())]

    # load the saved negative pairs for train + valid
    data_subfolder = data_dir + name + "_prepare_2/neg_pairs_test_hla_" + str(test_hla_i)

    hla_i_list = list(range(1, test_hla_i)) + list(range(test_hla_i + 1, limit_hla_i + 1))

    negative_train_valid = []
    # loop through all HLAs in train and valid
    for hla_i in hla_i_list:
        cur_fname = "neg_pairs_" + str(hla_i) + ".csv"
        cur_df = pd.read_csv(data_subfolder + "/" + cur_fname, header=0)
        cur_neg_pairs = [(tcr, hla) for tcr, hla in zip(cur_df.tcr.to_list(),
                                                        cur_df.hla.to_list())]
        negative_train_valid += cur_neg_pairs

    # compute the number of positive/negative for training
    np_train = round(len(positive_train_valid)*prop_train)
    nn_train = np_train * n_fold

    # shuffle pos and neg pairs
    shuffle(positive_train_valid)
    shuffle(negative_train_valid)
    # split positive pair list
    positive_train_ori = positive_train_valid[:np_train]
    positive_valid_ori = positive_train_valid[np_train:]
    # split negative pair list
    negative_train_ori = negative_train_valid[:nn_train]
    negative_valid_ori = negative_train_valid[nn_train:]
    # convert to dataframes
    df_train_pos = pd.DataFrame(positive_train_ori, columns=['tcr', 'hla_allele'])
    df_train_neg = pd.DataFrame(negative_train_ori, columns=['tcr', 'hla_allele'])
    df_valid_pos = pd.DataFrame(positive_valid_ori, columns=['tcr', 'hla_allele'])
    df_valid_neg = pd.DataFrame(negative_valid_ori, columns=['tcr', 'hla_allele'])
    # prepare folder for writing files out to
    model_data_dir = data_dir + name + "_raw_v"
    if not os.path.isdir(model_data_dir):
        os.mkdir(model_data_dir)
    # write out pos/neg train/valid
    df_train_pos.to_csv(
      model_data_dir + "/" + name + "_train_pos_" + str(test_hla_i) + ".csv", index=False)
    df_train_neg.to_csv(
      model_data_dir + "/" + name + "_train_neg_" + str(test_hla_i) + ".csv", index=False)
    df_valid_pos.to_csv(
      model_data_dir + "/" + name + "_valid_pos_" + str(test_hla_i) + ".csv", index=False)
    df_valid_neg.to_csv(
      model_data_dir + "/" + name + "_valid_neg_" + str(test_hla_i) + ".csv", index=False)


if __name__ == "__main__":
    args = parser.parse_args()
    print(args)
    get_data(args)
