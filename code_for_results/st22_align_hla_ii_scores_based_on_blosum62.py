# to prepare for computeing the distance matrix based on blosum62
# following the formula in NetMHCpan paper,
# first compute the alignment scores of the HLA pseudo sequences

import numpy as np
import pandas as pd

from collections import defaultdict
from collections import Counter

# load blosum62_X
blosum62_X = pd.read_csv("../data/for_encoders/blosum62_X.csv", header=0)
blosum62_X = blosum62_X.iloc[: , 1:]
blosum62_X.shape

# load HLA pseudo sequence
HLA_II_seq = pd.read_csv("../data/for_encoders/HLA_II_pseudo_45.csv", header=0)
HLA_II_seq.shape

# load HLA names in chosen order
# (following the same order as those used in predicted scores )

df_hla = pd.read_csv("../data/Liu_2019/liu_2019_hla_ii_format.csv", header=0)
df_hla.shape

aa_list = list(blosum62_X.columns)

blosum62_X_dict = defaultdict(int)
#blosum62_X.iloc[0, 5]

for i, aa in enumerate(aa_list):
    for j, bb in enumerate(aa_list):
        blosum62_X_dict[(aa, bb)] = blosum62_X.iloc[i, j]


# compute alignment score each each pair of pseudo sequences
def score_seqs(seq1, seq2):
    seq1_parse = list(seq1)
    seq2_parse = list(seq2)
    value = 0
    for _, (aa, bb) in enumerate(zip(seq1_parse, seq2_parse)):
        value += blosum62_X_dict[(aa, bb)]
    return value



# build a dict for hla and pseudo sequence
seq_dict = defaultdict(str)

for cur_hla, cur_seq in zip(HLA_II_seq.hla.tolist(), HLA_II_seq.seq.tolist()):
    seq_dict[cur_hla] = cur_seq

pseudo_seqs = [seq_dict[x] for x in df_hla.hla.tolist()]

score_lists = []
for i, seq1 in enumerate(pseudo_seqs):
    cur_list = []
    for j, seq2 in enumerate(pseudo_seqs):
        cur_list += [score_seqs(seq1, seq2)]
    score_lists += [cur_list]


df_score = pd.DataFrame(score_lists, columns = df_hla.hla.to_list())


df_score.to_csv("../results/st22_blosum62_X_seq_align_score_HLA_II.csv", index = False)
