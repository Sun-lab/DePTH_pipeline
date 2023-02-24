# this file creates a TCR appearance dictionary
# with the key being a TCR in the format of (CDR3, v_allele)
# and the value being a set of individuals that this TCR appears in
# from the data files of 666 individuals that DeWitt_2018 used

# only consider those with frame_type == 'In'

import os
import pickle
import numpy as np
import pandas as pd

from collections import defaultdict
from collections import Counter



directory = '../../data/emerson-2017-natgen/HIP_folder'
file_list = []
for entry in os.scandir(directory):
    file_list += [entry.path]

file_list.sort()

# create a dictionary to store where each (cdr3, v_allele) pair appears
indiv_dict = defaultdict(set)


# want to keep two kind of TCRs:
# part 1 with v_allele being not na
# part 2 with v_allele being na but v_gene not unresolved
#        - this part should be equivalent to the part where v_allele_ties
#          is not NaN
#        - in this case, take the first allele number in the ties

for i, file_name in enumerate(file_list):
    # load file
    cur = pd.read_csv(file_name, sep = "\t", header = 0, low_memory=False)
    # subset the data
    cur_in = cur[cur['frame_type'] == 'In']
    # for every TCR in the format of (CDR3, v_allele) with available or possible
    # options for v_allele,
    # add index i to the dictionary value set for this (CDR3, v_allele) tuple
    for aa, v_gene, v_allele, v_allele_ties in zip(cur_in.amino_acid.tolist(), \
                                                   cur_in.v_gene.tolist(), \
                                                   cur_in.v_allele.tolist(), \
                                                   cur_in.v_allele_ties.tolist()):
        # a trick making use of np.nan != np.nan
        if v_allele == v_allele:
            if int(v_allele) < 10:
                str_v_allele = '0' + str(int(v_allele))
            else:
                str_v_allele = str(int(v_allele))
            v_full_name = v_gene + '*' + str_v_allele
            indiv_dict[(aa, v_full_name)].add(i)
        else:
            if v_allele_ties == v_allele_ties:
                first_allele = v_allele_ties.split(',')[0]
                if int(first_allele) < 10:
                    str_v_allele = '0' + str(int(first_allele))
                else:
                    str_v_allele = str(int(first_allele))
                v_full_name = v_gene + '*' + str_v_allele
                indiv_dict[(aa, v_full_name)].add(i)


with open('../../data/intermediate_data/t5_in_allele_level_tcr_dict.pickle', 'wb') as f:
    pickle.dump(indiv_dict, f)
