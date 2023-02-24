This file summarizes what four jupyter notebooks:

t4_extend_pseudo_seq_HLA_I.ipynb

t4_extend_pseudo_seq_HLA_II_alpha.ipynb

t4_extend_pseudo_seq_HLA_II_beta.ipynb

t4_combine_HLA_II_AB.ipynb

do and find out.


We have pseudo sequences for HLA alleles on 34 positions for HLA-I alleles ("../../data/intermediate_data/HLA_I_v2_pseudo_sub.csv"), 15 positions for HLA-II alpha chain and 19 positions for HLA-II beta chain ("../../data/intermediate_data/pseudosequence_2016_all_X.dat").

We want to extend the positions to also cover those additional positions got from t2_check_additional_pos_contacts.log.ipynb. This is done in the first three separate jupyter notebooks. The fourth jupyter notebook t4_combine_HLA_II_AB.ipynb combines the extended pseudo sequences for HLA-II A/B alleles together to form a combined pseudo sequence for each HLA-II pair.

In each of the first three notebooks, we do the following steps:

First, we need to check whether we can get an unique pseudo sequence for each HLA on the additional positions.

If the unique sequences are verified in the first step, secondly, we will keep those positions with diversity in amino acids, and write out the pseudo sequences into files.

The resource data files for the full sequences:

"../../data/intermediate_data/HLA_TCR_contact/ClassI_prot.alfas"

"../../data/intermediate_data/HLA_TCR_contact/DRA_prot.alfas"

"../../data/intermediate_data/HLA_TCR_contact/DRB_prot.alfas"

"../../data/intermediate_data/HLA_TCR_contact/DPA1_prot.alfas"

"../../data/intermediate_data/HLA_TCR_contact/DPB1_prot.alfas"

"../../data/intermediate_data/HLA_TCR_contact/DQA1_prot.alfas"

"../../data/intermediate_data/HLA_TCR_contact/DQB1_prot.alfas"

The 34 positions (1-indexed) for HLA-I alleles from NetMHCPan are:

(7, 9, 24, 45, 59, 62, 63, 66, 67, 69, 70, 73, 74, 76, 77, 80, 81, 84, 95, 97, 99, 114, 116, 118, 143, 147, 150, 152, 156, 158, 159, 163, 167, 171)

The 15 positions (1-indexed, but needs other additional adjustment) for HLA-II A alleles from NetMHCIIpan-3.0 are:

(9, 11, 22, 24, 31, 52, 53, 58, 59, 61, 65, 66, 68, 72, 73)

The 19 positions (1-indexed, but needs other additional adjustment) for HLA-II B alleles from NetMHCIIpan-3.0 are:

(9, 11, 13, 26, 28, 30, 47, 57, 67, 70, 71, 74, 77, 78, 81, 85, 86, 89, 90)

The details for the additional adjustments can be found from the summary file for step70_expore_HLA-II.ipynb.

The additional position candidates (index ajusted) are:

17 additional positions for HLA-I: (4, 57, 64, 67, 71, 74, 122, 145, 148, 150, 153, 154, 156, 160, 161, 165, 169).

18 additional positions for HLA-II alpha chain: (4, 28, 35, 39, 45, 46, 47, 50, 51, 53, 56, 58, 60, 63, 65, 67, 71, 72).

9 additional positions for HLA-II beta chain: (49, 53, 54, 57, 59, 62, 66, 69, 75).


------------------------------------------------------------
Summary for t4_extend_pseudo_seqs_HLA_I.ipynb:
This code does four things:

(0) It verifies that for each of the 85 HLA-I alleles in v2 of DeWitt_2018 HLA matrix data, the pseudo sequences we can get on the 34 original positions based on full sequences in "../../data/intermediate_data/HLA_TCR_contact/ClassI_prot.alfas" are the same, if multiple corresponding full sequences exist. It also verifies that these pseudo sequences we reconstructed from full sequences match those in file "../../data/intermediate_data/HLA_I_v2_pseudo_sub.csv".

(1) It verifies that for each of the 85 HLA-I alleles, we can also find a unique pseudo sequence consisting of amino acids on the 17 positions.

(2) It finds out the 6 out of these 17 positions that have diversity in terms of amino acids among the 85 HLA-I alleles.

(3) It combines these 6 additional positions with the original 34 positions, and writes out a file of pseudo sequences based on the in total 40 positions.

Output file:

"../../data/for_encoders/step75_HLA_I_v2_pseudo_40.csv".

--------------------------------------------------------------

Summary for t4_extend_pseudo_seqs_HLA_II_alpha.ipynb:

The summary for this part contains lots of details due to the adjustments needed to do for alignments. So we start from stating the findings on how to adjust the positions to get the pseudo sequences.

(0) In general, the 15 positions (9, 11, 22, 24, 31, 52, 53, 58, 59, 61, 65, 66, 68, 72, 73) given in NetMHCIIpan-3.0 for HLA-II alpha chains need to be subtracted by 4 in order to match the info given in full sequences files(DRA_prot.alfas, DPA_prot.alfas, DQA_prot.alfas).

(1) On top of that, the first one (pos 9) of these 15 positions needs to be subtracted by an additional 1 (one guess is it is related to a deletion in 0-indexed position 5, for example, as in the unique full sequence 'HVIIQ-AEFYLNPDQSGEFMFDFDGDEIFHVDMAKKETVWRLEEFGRFASFEAQGALANIAVDKANLEIMTKRSN' for DRA alleles, and the alignment in NetMHCIIpan-3.0 might have ignored this deletion and caused the relative distance between the first postion and the second to be shorter by 1).

(2) Besides these, there is an additional modification needed to be done if the alpha chains falls into a set of DQAs: ('DQA1\*05:04','DQA1\*05:07', 'DQA1\*05:01', 'DQA1\*05:06', 'DQA1\*05:03', 'DQA1\*05:02', 'DQA1\*05:08', 'DQA1\*05:05', 'DQA1\*05:11', 'DQA1\*05:10', 'DQA1\*05:09', 'DQA1\*04:02', 'DQA1\*04:04', 'DQA1\*06:01', 'DQA1\*06:02', 'DQA1\*02:01', 'DQA1\*04:01'). Most of these DQAs are got from Fig.2 of NetMHCIIpan-3.0 paper, except that 'DQA1\*04:01' was not mentioned in the paper but found out to have the same behavior in terms of adjusting positions for pseudo sequence matching, and so it is added into this list. Based on the paper, the sequences of these HLA-II alpha chains have a deletion on position 53 (under indexing system in NetMHCIIpan-3.0 paper). However, the full sequences from DQA_prot.alfas show no deletion here but a deletion on position 48 (under indexing system in NetMHCIIpan-3.0 paper).
For example, DQA1_05_04 has sequence (starting from position 35)

DLGRKETVWCLPVLRQFR-FDPQFALTNIAVLKHNLN...
v.s.

DLGRKETVWCLPV-LRQFRFDPQFALTNIAVLKHNLN...

given by full sequene info as in DQA_prot.alfas.

Thus, the deletion at postion 48 is not considered by the alignment of NetMHCIIpan-3.0 as a deletion, but the one at position 53 is. So in order to construct pseudo sequence in the same way as in pseudosequence_2016_all_X.dat, for these special DQAs, one additional adjustment we need to make is that, to amino acids for positions 48, 49, 50, 51, 52, we need to take those on positions 49, 50, 51, 52, 53 (under indexing system of NetMHCIIpan-3.0) from the full sequences in DQA_prot.alfas instead, and if we need to get amino acid on position 53, write an "X" instead. For other DQAs that do not fall into this special set, DPAs and the DRA, this adjustment should not be done.

(3) If one HLA-II A allele has multiple corresponding pseudo sequences, then ignore those sequences starting with "X" and only keep the others.


Besides finding out these adjustments for alignments, this code does the following things:

(0) There are 17 HLA-II alpha chains (here we treat all DRAs as one, since they all have the same full amino acid sequences) involved in the HLA-II pairs in v2 of DeWitt_2018 HLA matrix data.

(1) For each of these 17 HLA-II alpha chains, based on full sequence info (DRA_prot.alfas, DPA1_prot.alfas, DQA1_prot.alfas) and the original 15 positions from NetMHCIIpan-3.0, after adjustments to the position indexes, we can reconstruct a unique pseudo sequence, and the code verifies that the reconstructed pseudo sequence matches that from pseudosequence_2016_all_X.dat.

(2) Out of the 18 additional positions, the first one is position 8. This position is ignored as it is the same as the first one (9) in the original 15 positions, after that 9 is shifted by 1 to the left and becomes 8. For the remaining 17 positions, we make the adjustments mentioned above and get the amino acids of the 17 HLA-II A alleles on these positions. Out of these 17 positions, we find out 7 with diversiy in terms of amino acids.

(3) The 7 additional positions are combined with the original 15 positions into a set of 22 positions. Under the adjustments mentioned above, we construct pseudo sequences (verified to be unique for each allele) for each HLA-II A allele, and write the dictionary with allele key and pseudo sequence as value out into a csv file: "../../data/intermediate_data/t4_HLA_II_v2_alpha_pseudo_22_dict.csv".

--------------------------------------------------------------

Summary for t4_extend_pseudo_seqs_HLA_II_beta.ipynb:

The summary for this part contains details due to the adjustments needed to do for alignments. We start from stating the findings on how to adjust the positions to get the pseudo sequences.

(0) In general, the 19 positions (9, 11, 13, 26, 28, 30, 47, 57, 67, 70, 71, 74, 77, 78, 81, 85, 86, 89, 90) given in NetMHCIIpan-3.0 for HLA-II alpha chains need to be subtracted by 7 in order to match the info given in full sequences files(DRB_prot.alfas, DPB1_prot.alfas, DQB1_prot.alfas).

(1) Besides that, there is an additional modification needed to be done if the B allele is either 'DPB1\*10:01' or 'DPB1\*17:01'. In this case, the first one (pos 9) of these 19 positions needs to be subtracted by an additional 2 (for this one I do not have a reasonable guess for the cause).
One guess would be that it is caused by deletion. However,

The full sequence of 'DRB1*01:01' is

'FLWQLKFECHFFNGTERVRLLERCIYNQEESVRFDSDVGEYRAVTELGRPDAEYWNSQKDLLEQRRAAVDTYCRHNYGVGESFTV'

The full sequence of 'DPB1*01:01' is

'YVYQGRQECYAFN--GTQRFLERYIYNREEYARFDSDVGEFRAVTELGRPAAEYWNSQKDILEEKRAVPDRVCRHNYELDEAVTL'

The full sequence of 'DPB1*02:01' is

'YLFQGRQECYAFN--GTQRFLERYIYNREEFVRFDSDVGEFRAVTELGRPDEEYWNSQKDILEEERAVPDRMCRHNYELGGPMTL'

The full sequence of 'DPB1*10:01' is

'YVHQLRQECYAFN--GTQRFLERYIYNREEFVRFDSDVGEFRAVTELGRPDEEYWNSQKDILEEERAVPDRVCRHNYELDEAVTL'

The full sequence of 'DPB1*17:01' is

'YVHQLRQECYAFN--GTQRFLERYIYNREEFVRFDSDVGEFRAVTELGRPDEDYWNSQKDILEEERAVPDRMCRHNYELDEAVTL'

Thus all 'DPB1\*01:01', 'DPB1\*02:01', 'DPB1\*10:01' and 'DPB1\*17:01' have deletions, but neither 'DPB1\*01:01' nor 'DPB1\*02:01' needs this additional adjustment. For 'DPB1\*02:01', the first aa is 'Y' and is different from the third that is 'F', so it is not a coincident of having the same amino acid on two positions before and after the shift (pos 2 and 0 under 0-indexing).

(2) If one HLA-II A allele has multiple corresponding pseudo sequences, then ignore those sequences containing with "X" and only keep the others.


Besides these adjustments for alignments, this code does the following things:

(0) There are 62 HLA-II B alleles involved in the HLA-II pairs in v2 of DeWitt_2018 HLA matrix data.

(1) For each of these 62 HLA-II B alleles, based on full sequence info (DRB_prot.alfas, DPB1_prot.alfas, DQB1_prot.alfas) and the original 19 positions from NetMHCIIpan-3.0, after adjustments to the position indexes, we can reconstruct a unique pseudo sequence, and the code verifies that the reconstructed pseudo sequence matches that from pseudosequence_2016_all_X.dat.

(2) Based on amino acids from the full sequences on the 9 additional positions after shifting by 7, we verify that we can construct a unique pseudo sequence for each of the 62 HLA-II B alleles. Out of these 9 positions, we find out 4 with diversiy in terms of amino acids.

(3) The 4 additional positions are combined with the original 19 positions into a set of 23 positions. Under the adjustments mentioned above, we construct pseudo sequences (verified to be unique for each allele) for each HLA-II B allele, and write the dictionary with allele key and pseudo sequence as value out into a csv file: "../../data/intermediate_data/t4_HLA_II_v2_beta_pseudo_23_dict.csv".


--------------------------------------------------------------

Summary for t4_combine_HLA_II_AB.ipynb

For each of 135 HLA-II A/B allele pairs from HLA v2 matrix data of DeWitt_2018 (the original data has 130 items including 5 haplotypes which are each separated into 2 pairs, resulting in 135 pairs in total), this file gets the corresponding A & B alleles and combines the corresponding two extended pseudo sequences (A len 22, B len 23) together to form one pseudo sequence of length 45.

The input files are:

"../../data/intermediate_data/DeWitt_2018/HLA_v2_features_row_names.txt"

"../../data/intermediate_data/t4_HLA_II_v2_alpha_pseudo_22_dict.csv"

"../../data/intermediate_data/t4_HLA_II_v2_beta_pseudo_23_dict.csv"

The output file is:

"../../data/for_encoders/HLA_II_pseudo_45.csv".
