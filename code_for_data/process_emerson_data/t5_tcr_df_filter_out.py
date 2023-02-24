# (1) convert allele level tcr dict into data frame
# after filtering out the tcr tuple that showed up in
#   "../data/elife-38358-fig11-data1-v2.tds
# (2) create another file containing only public tcrs

import os
import re
import pickle
import numpy as np
import pandas as pd

from collections import defaultdict
from collections import Counter


data_dir = '../../data/intermediate_data/'

allele_tcr_dict = pickle.load(open(data_dir+'step40_in_allele_level_tcr_dict.pickle', 'rb'))
#len(allele_tcr_dict)
#69122667


# filtered out tcr info
# note that the format of v_family from "elife-38358-fig11-data1-v2.tds"
#   is different from that in allele_tcr_dict
filter_tcr_info = pd.read_csv(data_dir+"DeWitt_2018/elife-38358-fig11-data1-v2.tds",
                              sep = "\t", header = 0)
filter_tcr_info.nunique()
# V-family             25
# CDR3-amino-acids    205
# CDR3-nucleotides    276
# Occurrence-count     77
# P_gen               531
filtered_key_list = [(v_family, aa) for v_family, aa in \
          zip(filter_tcr_info['V-family'], filter_tcr_info['CDR3-amino-acids'])]
filtered_key_set = set(filtered_key_list)
# no strange v_famiy show up in filter_tcr_info
Counter(filter_tcr_info['V-family'])
set(filter_tcr_info['CDR3-amino-acids'])

# one column for v_family
# one column for aa
# one column counting the # of individuals
# one column for a string of individuals
v_allele_col = []
aa_col = []
n_ind = []
string_ind = []

public_v_allele_col = []
public_aa_col = []
public_n_ind = []
public_string_ind = []

full_allele_key_list = list(allele_tcr_dict.keys())

for allele_key in full_allele_key_list:
    aa = allele_key[0]
    v_allele = allele_key[1]
    v_family = re.split('-|\*', v_allele)[0][4:]
    if (v_family, aa) not in filtered_key_set:
        v_allele_col += [v_allele]
        aa_col += [aa]
        value = list(allele_tcr_dict[allele_key])
        n_ind += [len(value)]
        value.sort()
        string_ind += [','.join([str(ind) for ind in value])]
        if len(value) > 1:
            public_v_allele_col += [v_allele]
            public_aa_col += [aa]
            public_n_ind += [len(value)]
            public_string_ind += [','.join([str(ind) for ind in value])]

df_allele_level_tcr = pd.DataFrame(\
                           list(zip(v_allele_col, aa_col, n_ind, string_ind)),
                           columns = ['v_allele', 'amino_acids', 'n_ind', 'individuals'])
#df_allele_level_tcr.shape
# (69121789, 4)
df_allele_level_tcr.to_csv(data_dir+"t5_allele_level_tcr_filtered_wrt_vf_and_aa.csv",
                           index = False)


df_public_allele_level_tcr = pd.DataFrame(\
                                 list(zip(public_v_allele_col, public_aa_col, public_n_ind, public_string_ind)),
                                 columns = ['v_allele', 'amino_acids', 'n_ind', 'individuals'])
#df_public_allele_level_tcr.shape
# (8739207, 4)
df_public_allele_level_tcr.to_csv(\
    data_dir+"t5_public_allele_level_tcr_filtered_wrt_vf_and_aa.csv", index = False)







#
