
Construct data sets for HLA-I leave1out experiments

The purpose is to make sure that only the chosen HLA-I allele is allowed in test pairs and no TCR is allowed to overlap between training + validation and test pairs. Subject to this constraint, want to maximize the number of positive pairs that can be used.

[sample_neg_tcrs_leave1out_match_test.py]

For a given HLA-I allele reserved for test purpose, use all positive pairs involving this HLA-I allele as positive test pairs, filter the positive pairs left for training and validation by excluding pairs involving TCRs that appear in the test positive pairs, sample TCRs for negative pairs in the same frequency matching manner as for HLA_I_all_match, except that any TCRs involved in the positive training and validation pairs are also excluded. The index at the end of output files like HLA_I_leave1out_match_{number}.csv corresponds to the index (1-indexed) of the reserverd HLA-I alleles among all HLA-I alleles in the table HLA_v2_features_reformat.csv, with 1 for HLA-B\*08:01, 11 for HLA-B\*07:02, and 23 for HLA-C\*07:01.

[sample_neg_tcrs_leave1out_match_train_valid.py]

Under each setting of reserving one HLA-I allele for test data, for each of the other 84 HLA-I alleles, sample TCRs to form negative pairs with, in order to prepare for negative pairs in training and validation data sets. The sampling is still done in the manner of frequency matching, but this time no TCR in either test positive or test negative pairs is allowed to appear in the sampled TCRs.

[verify_sufficient_negative_TCRs.ipynb]

verify there are enough unique negative TCRs to sample from based on one output from the step above. The answer is yes for all three settings.

[generate_data.py]

Under each setting of reserving one HLA-I allele for test data, form the sets of all positive pairs and the set of sampled negative pairs, and split into 75% training and 25% validation (note that the relative proportion between training/validation and test is not the same as before, since all positive pairs involving the reserved HLA allele are included in test, and training/validation positive pairs can only use other positive pairs involving the other HLA-I alleles after some filtering criteria on TCRs). run by .sh file [generate_data_leave1out.sh].

[hla_i_all_leave1out_match_translate.ipynb]

translate v gene format to the format as in file ../../../data/for_encoders/combo_xcr.tsv and reorganize files.
