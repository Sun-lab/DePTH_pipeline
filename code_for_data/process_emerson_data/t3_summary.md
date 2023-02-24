This file is a summary for t3_expore_HLA-II.ipynb.

t3_expore_HLA-II.ipynb works on verifying the pseudo sequence of HLA-II pairs (each pair consists of an alpha chain and a beta chain) from HLA_v2_features.txt of DeWitt_2018.

Data files involved:

"../../data/intermediate_data/DeWitt_2018/pseudosequence_2016_all_X.dat" (the pseudo sequence data of HLA-II pairs)

"../../data/intermediate_data/DeWitt_2018/HLA_features_row_names.txt" (name of HLA alleles v1 from DeWitt_2018)

"../../data/intermediate_data/DeWitt_2018/HLA_v2_features_row_names.txt" (name of HLA v2 from DeWitt_2018)

(The six files below give the sequences for DRA, DRB, DPA, DPB, DQA and DQB alleles.)

"../../data/intermediate_data/HLA_TCR_contact/DRA_prot.alfas"
"../../data/intermediate_data/HLA_TCR_contact/DRB_prot.alfas"
"../../data/intermediate_data/HLA_TCR_contact/DPA_prot.alfas"
"../../data/intermediate_data/HLA_TCR_contact/DPB_prot.alfas"
"../../data/intermediate_data/HLA_TCR_contact/DQA_prot.alfas"
"../../data/intermediate_data/HLA_TCR_contact/DQB_prot.alfas"

(The two lines below are important positions on HLA from the paper:
NetMHCIIpan-3.0, a common pan-specific MHC class II prediction method including all three human MHC class II isotypes, HLA-DR, HLA-DP and HLA-DQ.)

pos_alpha = [9, 11, 22, 24, 31, 52, 53, 58, 59, 61, 65, 66, 68, 72, 73]
pos_beta = [9, 11, 13, 26, 28, 30, 47, 57, 67, 70, 71, 74, 77, 78, 81, 85, 86, 89, 90]



What this jupyter notebook does and finds out in more details:


HLA_features v1 contains 173 HLA alleles while HLA_v2_features contains 215 (including 85 HLA-I alleles and 130 items for HLA-II, which contain: DRB, DPA & DPB appearing in pairs, DQA & DQB appearing in pairs, and additional 5 DR-DQ haplotypes).

The components of all these 130 items related to HLA-II in v2 are verified to exist in v1 in the format of single HLA-II alleles.

To verify the pseudo sequence alignment, we break the 5 DR-DQ haplotypes into 5 DRB and 5 DQ pairs and get a set of 135 unique HLA-II pairs (each DRB actually corresponds to a (DRA, DRB) pair, since the sequence for DRA alleles is unique). This set contains 25 DP pairs, 72 DQ pairs and 39 DRB1 pairs.

"../../data/intermediate_data/HLA_TCR_contact/DRA_prot.alfas" contains sequences for 7 DRA alleles (to the level of two colons or three colons), but the 7 sequences are all the same ('HVIIQ-AEFYLNPDQSGEFMFDFDGDEIFHVDMAKKETVWRLEEFGRFASFEAQGALANIAVDKANLEIMTKRSN'). If we align it with the sequences shown in figure 2 on page 4 of the NetMHCIIpan-3.0 paper, the position 31 (0-indexed) of the sequence from DRA_prot.alfas corresponds to the position 35(1-indexed) in that from NetMHCIIpan-3.0.

"../../data/intermediate_data/HLA_TCR_contact/DRB_prot.alfas" contains sequences for 2238 DRB alleles (not unique if we cut them off after the second colon). If cut off after the second colon, with the resulting short name as the key and corresponding set of sequences as the value, we get a dictionary of size 1709, with 1697 keys having 1 corresponding sequences, 11 having two unique sequences and 1 having 3.

"../../data/intermediate_data/HLA_TCR_contact/DPA1_prot.alfas" contains sequences for 53 DQA alleles (not unique if we cut them off after the second colon). If cut off after the second colon, with the resulting short name as the key and corresponding set of sequences as the value, we get a dictionary of size 24.

"../../data/intermediate_data/HLA_TCR_contact/DPB1_prot.alfas" contains sequences for 873 DQA alleles (not unique if we cut them off after the second colon). If cut off after the second colon, with the resulting short name as the key and corresponding set of sequences as the value, we get a dictionary of size 641.

"../../data/intermediate_data/HLA_TCR_contact/DQA1_prot.alfas" contains sequences for 91 DQA alleles (not unique if we cut them off after the second colon). If cut off after the second colon, with the resulting short name as the key and corresponding set of sequences as the value, we get a dictionary of size 35.

"../../data/intermediate_data/HLA_TCR_contact/DQB1_prot.alfas" contains sequences for 1045 DQA alleles (not unique if we cut them off after the second colon). If cut off after the second colon, with the resulting short name as the key and corresponding set of sequences as the value, we get a dictionary of size 748.


When looking into one example per allele type (DRA, DRB, DPA, DPB, DQA and DQB) to find the matching positions for verifying consistency, here are the findings:

DRA's unique sequence matches that in the pseudo sequences if we shift the pos_alpha indexes all by subtracting 4 (this is from 1-indexed to 0-indexed, so it would corresponds to shifting by 3 to the left if still under 1-indexed setting), and replace the first resulting aa by the one before it (my guess is the difference might be due to some deletion).

DRB sequence dictionary has 1 unique sequence for 1697 alleles (cutoff after two colons), 2 sequences for 11 alleles and 3 for 1 allele. The example pseudo sequence was matched after we subtract the pos_beta indexes all by 7 (should be moving 6 positions to the left if still under 1-indexed setting).

DPA sequence dictionary has 1 unique sequence for 21 alleles (cutoff after two colons) and 2 for 3. The sequence match is found through subtracting the indexes in pos_alpha by 4, and replacing the first resulting aa by the one before it, but will need to take the one without "X" as aa, if there are more than one unique sequences.

DPB sequence dictionary has 1 unique sequence for 640 alleles (cutoff after two colons) and 2 for 1 allele. The example pseudo sequence was matched after we subtract the pos_beta indexes all by 7.

DQA sequence dictionary has 1 unique sequence for 34 alleles (cutoff after two colons) and 2 for 1. The sequence match is found through subtracting the indexes in pos_alpha by 4, and replacing the first resulting aa by the one before it, but might need to take the one without "X" as aa, if there are more than one unique sequences.

DQB sequence dictionary has 1 unique sequence for 743 alleles (cutoff after two colons) and 2 for 5 alleles. The example pseudo sequence was matched after we subtract the pos_beta indexes all by 7. But might need to take the one without "X" as aa, if there are more than one unique sequences.

Finally move to verify the consistency.

First try on directly matching using following the principles below has only 74.8% match:
(0) DPA, DQA are processed by subtracting 4 from the indexes from pos_alpha, with the first aa replaced by the one before it;
(1) If multiple sequences are found, take the subsequences according to the indexes and abandon those containing "X". If there are more than one unique subsequences left or no left, raise an error flag.

It seems that out of those that match, the DQA parts all correspond to those with deletion in figure 2 on page 4 of NetMHC_II_pan-3.0 paper. Out of those that do not match, the DQA parts all correspond to those without deletion.

Looking into one specific allele DQA1\*05:01, the sequence from DQA_prot.alfas is:
'HVASYGVNLYQSYGPSGQYTHEFDGDEQFYVDLGRKETVWCLPV-LRQFRFDPQFALTNIAVLKHNLNSLIKRSN'
and the one from NetMHC_II_pan-3.0 paper is:
'DLGRKETVWCLPVLRQFR-FDPQFALTNIAVLKHNLNSLIKRSNSTAATN----------'
My understanding is that it is due to different ways of aligning the sequences. The deletion on position 53 from NetMHC_II_pan-3.0 paper is not considered as a deletion in DQA_prot.alfas, but a delection is considered by DQA_prot.alfas before position 48 (using the index from NetMHC_II_pan-3.0). Thus due to that pos_alpha does not contain other positions in the area, the reconstructed pseudo sequences and those from "../../data/intermediate_data/pseudosequence_2016_all_X.dat" still match on other positions except on 52 & 53 (1-indexed under the setting of NetMHC_II_pan-3.0).

For the 4 DP pairs which do not have the pseudo sequences matching the reconstructed ones, the trouble are on DPBs, specifically, DPB1\*10:01 and DPB1\*17:01. DPB1\*10:01 is supposed to be 'YLQFEYFDIEEVRVHLDVT' based on DPB_prot.aflas, but turns out to be 'YLQFEYFDIEEVRVHLDVT', with difference in the first aa. DPB1\*17:01 is supposed to be 'HLQFEYFDIEEVRMHLDVT', but turns out to be 'YLQFEYFDIEEVRMHLDVT', again with difference in the first aa. The reason is not known yet.

The way to fix the verification procedure for now:

(0) If DQA falls into the list extra_modify_DQAs, need to modify positions 52 & 53 (1-indexed under NetMHC_II_pan-3.0);
(1) If DPB falls into the list extra_modify_DPBs, need to manually modify the first position to 'Y'.
(2) Change the getting rid of all pseudo sequences containing 'X' condition inside function get_a_half_modify to only getting rid of those starting with 'X'.

This time all 135 pairs each has pseudo sequence from "../../data/intermediate_data/pseudosequence_2016_all_X.dat" matching that reconstructed from prot.alfas files.
