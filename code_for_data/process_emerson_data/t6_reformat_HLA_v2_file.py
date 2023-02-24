# this file prepares HLA_v2 input file for computing p-values in
# later steps

# Input file:
#   ../../data/intermediate_data/DeWitt_2018/HLA_v2_features.txt
# Output file:
#   ../../data/intermediate_data/HLA_v2_features_reformat.csv

import os
import re
import pickle
import numpy as np
import pandas as pd

from collections import defaultdict
from collections import Counter


# HLA information

name_list = []
n_pos_list = []
ind_pos_list = []
n_neg_list = []
ind_neg_list = []

f = open("../../data/intermediate_data/DeWitt_2018/HLA_v2_features.txt", "r")

for line in f:
    line_list = re.split(' |\n', line)
    name = line_list[1]
    n_pos = int(line_list[3])
    ind_pos = line_list[5:(5+n_pos)]
    n_neg = int(line_list[6+n_pos])
    ind_neg = line_list[(8+n_pos):(-1)]
    # add them to list
    name_list += [name]
    n_pos_list += [n_pos]
    ind_pos_list += [ind_pos]
    n_neg_list += [n_neg]
    ind_neg_list += [ind_neg]

n_total = [n_pos + n_neg for n_pos, n_neg in zip(n_pos_list, n_neg_list)]
max(n_total)
#630
min(n_total)
#466

ind_pos_string_list = [','.join(inds) for inds in ind_pos_list]
ind_neg_string_list = [','.join(inds) for inds in ind_neg_list]
# write data frame out for the purpose of running R code
HLA_v2_features_reformat = \
       pd.DataFrame(list(zip(name_list, n_pos_list, ind_pos_string_list, \
                             n_neg_list, ind_neg_string_list)), \
                             columns = ['hla', 'n_pos', 'ind_pos', \
                                        'n_neg', 'ind_neg'])
HLA_v2_features_reformat.to_csv("../../data/intermediate_data/HLA_v2_features_reformat.csv", \
                                 index = False)


#
