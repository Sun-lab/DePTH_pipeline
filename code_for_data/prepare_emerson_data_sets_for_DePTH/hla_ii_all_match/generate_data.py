#!/usr/bin/env python3

import os
import sys
import pandas as pd

import random
from random import shuffle

from collections import Counter
from collections import defaultdict

import argparse

# the following line is to deal with the case that self-defined module
# in the submission directory are not found when submitted through sbatch, as in
# https://stackoverflow.com/questions/46718499/
# how-does-one-make-sure-that-the-python-submission-script-in-slurm-is-in-the-loca/
# 46724189?noredirect=1#comment80425963_46724189
sys.path.append(os.getcwd())


parser = argparse.ArgumentParser(description='Generate (TCR, HLA) pairs for building the model')
parser.add_argument('--hla_class', default='HLA_I', help='the class of HLA')
parser.add_argument('--label', default='all_match',
                    help='part of the label of the model to tell different settings apart')
parser.add_argument('--prop_train', type=float, default='0.6',
                    help='proportion of positive pairs to use for training')
parser.add_argument('--prop_valid', type=float, default='0.2',
                    help='proportion of positive pairs to use for validation')
parser.add_argument('--n_fold', type=int, default=5,
                    help='number of negative pairs divided by that of positive pairs')
parser.add_argument('--rseed', type=int, default=1627, help='random seed for random')


# generate pos/neg pairs for train/valid/test

def get_data(hla_class="HLA_I", label="all_random",
             prop_train=0.6, prop_valid=0.2, n_fold=5, rseed=1627):

    hla_class = args.hla_class
    label = args.label
    prop_train = args.prop_train
    prop_valid = args.prop_valid
    n_fold = args.n_fold
    rseed = args.rseed

    n_fold = int(n_fold)
    rseed = int(rseed)

    data_dir = "../../../data/intermediate_data/"
    # join hla class and part of the label to form a label for telling data resources apart
    name = hla_class + "_" + label
    # load TCR information
    TCR_rmat = pd.read_csv(data_dir + "public_allele_level_tcr_name.csv", header=0)
    TCR_name = [v + ',' + aa for v, aa in zip(TCR_rmat.v_allele.tolist(), TCR_rmat.amino_acids.tolist())]
    # move on to select the positive and negative pairs
    # load the association matrix with specified hla class
    if hla_class == "HLA_I":
        asso_file = data_dir + "t12_HLA_I_associated_TCR_v_alleles.csv"
    elif hla_class == "HLA_II":
        asso_file = data_dir + "t12_HLA_II_no_haplotype_associated_TCR_v_alleles.csv"
    else:
        print("hla_class input wrong. should be either HLA_I or HLA_II")
        return

    df_asso = pd.read_csv(asso_file, sep=",", header=0)
    # use dictionary to keep the associated pairs
    asso_dict = defaultdict(lambda: defaultdict(int))
    for tcr, hla in zip(df_asso.tcr, df_asso.hla_allele):
        asso_dict[hla][tcr] = 1
    # set random seed
    random.seed(rseed)
    # generate the pair list for positive and negative
    # positive:
    positive_ori = [(tcr, hla) for tcr, hla in zip(df_asso.tcr.to_list(), df_asso.hla_allele.to_list())]
    # compute the number of positives for train, valid and test
    np_train = round(len(positive_ori)*prop_train)
    np_valid = round(len(positive_ori)*prop_valid)
    np_test  = len(positive_ori) - (np_train + np_valid)
    # compute the number of negatives for train, valid and test
    nn_train = np_train * n_fold
    nn_valid = np_valid * n_fold
    nn_test = np_test * n_fold
    # different options for negative pairs
    if label == "all_random":
        # negative:
        shuffle(TCR_name)
        # count how many positive pairs each HLA appears in
        positive_hla_list = [item[1] for item in positive_ori]
        counter_pos_hla = Counter(positive_hla_list)
        positive_hla_unique = list(counter_pos_hla.keys())
        # construct the negative pairs that we really use
        negative_asso_dict = defaultdict(lambda: defaultdict(int))
        i = 0
        # negative test
        negative_ori = []
        for hla in positive_hla_unique:
            cnt_add = 0
            while cnt_add < (counter_pos_hla[hla]*n_fold):
                i = i % len(TCR_name)
                tcr = TCR_name[i]
                if tcr not in asso_dict[hla]:
                    if tcr not in negative_asso_dict[hla]:
                        negative_ori.append((tcr, hla))
                        negative_asso_dict[hla][tcr] = 1
                        cnt_add += 1
                i += 1

    elif label == "all_match":
        # negative are sampled to match the tcr frequency of those in positive
        negative_ori = []
        # add all negative pairs involving each HLA
        n_hlas = len(set(df_asso.hla_allele))
        for hla_i in range(1, n_hlas + 1):
            input_fname = "neg_pairs_" + str(hla_i) + ".csv"
            df_neg = pd.read_csv(data_dir + hla_class + "_all_match_prepare_2/"+input_fname, header = 0)
            cur_negs_ori = [(tcr, hla) for tcr, hla in zip(df_neg.tcr, df_neg.hla)]
            negative_ori += cur_negs_ori
    else:
        print("negative tcr selection method not coded for yet")
        return
    # shuffle pos and neg pairs
    shuffle(positive_ori)
    shuffle(negative_ori)
    # split positive pair list
    positive_test_ori = positive_ori[:np_test]
    positive_valid_ori = positive_ori[np_test:(np_test+np_valid)]
    positive_train_ori = positive_ori[(np_test+np_valid):]
    # split negative pair list
    negative_test_ori = negative_ori[:nn_test]
    negative_valid_ori = negative_ori[nn_test:(nn_test+nn_valid)]
    negative_train_ori = negative_ori[(nn_test+nn_valid):]
    # convert to dataframes
    df_train_pos = pd.DataFrame(positive_train_ori, columns=['tcr', 'hla_allele'])
    df_train_neg = pd.DataFrame(negative_train_ori, columns=['tcr', 'hla_allele'])
    df_valid_pos = pd.DataFrame(positive_valid_ori, columns=['tcr', 'hla_allele'])
    df_valid_neg = pd.DataFrame(negative_valid_ori, columns=['tcr', 'hla_allele'])
    df_test_pos  = pd.DataFrame(positive_test_ori,  columns=['tcr', 'hla_allele'])
    df_test_neg  = pd.DataFrame(negative_test_ori,  columns=['tcr', 'hla_allele'])
    # prepare folder for writing files out to
    output_dir = data_dir + name + "_raw_v"
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    # write out pos/neg train/valid/test
    df_train_pos.to_csv(output_dir + "/" + name + "_train_pos.csv", index=False)
    df_train_neg.to_csv(output_dir + "/" + name + "_train_neg.csv", index=False)
    df_valid_pos.to_csv(output_dir + "/" + name + "_valid_pos.csv", index=False)
    df_valid_neg.to_csv(output_dir + "/" + name + "_valid_neg.csv", index=False)
    df_test_pos.to_csv(output_dir  + "/" + name + "_test_pos.csv", index=False)
    df_test_neg.to_csv(output_dir  + "/" + name + "_test_neg.csv", index=False)


if __name__ == "__main__":
    args = parser.parse_args()
    print(args)
    get_data(args)
