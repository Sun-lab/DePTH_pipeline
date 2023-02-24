


The dataset that we use here was initially reported by

  Emerson, R. O., DeWitt, W. S., et al. (2017). Immunosequencing identifies signatures of cytomegalovirus exposure history and HLA-mediated effects on the T cell repertoire. Nature genetics, 49(5), 659-665.

  It was later amended (with additional HLA and covariate information) and re-analyzed by

  DeWitt III, W. S., Smith, A., Schoch, G., Hansen, J. A., Matsen IV, F. A., & Bradley, P. (2018). Human T cell receptor occurrence patterns encode immune history, genetic background, and receptor specificity. Elife, 7, e38358.

  We used the data reported by DeWitt (2018), and the file HLA_v2_features.txt is from pubtcrs_data_v1.tgz which was downloaded from Zenodo database (doi:10.5281/zenodo.1248193) on 10/02/2020.



# Prepare Emerson data HLA pseudo sequences

## Data exploration

[t1_extract_HLA_I_v2_pseudo.ipynb]

extracts the short pseudo sequences from MHC_pseudo.txt for 85 HLA-I alleles from DeWitt 2018.

HLA_v2_features.txt contains 215 (including HLA-I, HLA-II, and 5 DR-DQ haplotypes), among which 85 are HLA-I alleles.

MHC_pseudo.txt has 5\,327 rows which has 5\,315 unique names by the format (the 12 replicates are in Mamu and SLA, having nothing to do with HLA). These contains 4\,396 unique HLA alleles names for the format. All the 85 HLA-I alleles from HLA_v2_features.txt are included here. From MHC_pseudo.txt, we take out the subset corresponding to these 85 HLA-I alleles and save as file "HLA_I_v2_pseudo_sub.csv". (As a sanity check, for these 85 alleles, the rows in the original MHC_pseudo.txt with name formats HLA-\*\*\*:\*\* and HLA-\*\*\*\*\* actually have the same content in seq column, so these rows are duplicates and probably simply due to merging different formats)

## Explore possible additional contact positions

The jupyter notebook contains many details:

[t2_check_additional_pos_contacts.log.ipynb]

But a summary of what this notebook does and the findings is put in:

[t2_summary.md]


## Verify alignment info of HLA-II pairs in HLA data v2 from DeWitt_2018

The jupyter notebook doing this is long and contains many details:

[t3_expore_HLA-II.ipynb]

But a summary of what this notebook does and the findings is put in:

[t3_summary.md]


## Extend HLA contact positions

The first three jupyter notebook contains many details:

[t4_extend_pseudo_seqs_HLA_I.ipynb]

[t4_extend_pseudo_seqs_HLA_II_alpha.ipynb]

[t4_extend_pseudo_seqs_HLA_II_beta.ipynb]

[t4_combine_HLA_II_AB.ipynb]

But a summary of what these four notebooks do and the findings is put in:

[t4_summary.md]



# Prepare Emerson data associated TCR-HLA pairs.


## Construct data file of v_allele-level tcr existence among 666 individual files

The Emerson data used are downloaded from https://clients.adaptivebiotech.com/pub/Emerson-2017-NatGen. The repertoire data files are not included here since they are too big.

**Note:** the repertoire files used in this analysis were downloaded in year 2021, while we found that as of April 28, 2022, the format of the repertoire files in the link changed. The files were renamed and the columns are different from those before, although the information (after certain processing) can be matched to that in the old version.

[t5_organzie_repertoire_files.py]

organize Emerson repertoire files to put the HIP batch (666) in one folder and the Keck (120) batch in another.

[t5_make_allele_level_tcr_dict_from_666.py]

Construct a dictionary for allele level TCR appearance based on HIP batch of 666 subjects. The key of the dictionary is a TCR in the format of (CDR3, v_allele). The corresponding value is a set of the indexes of all the individuals that this TCR appears in. There are 69,122,667 unique TCRs. The output file is too big for github and is not included here.

[t5_tcr_df_filter_out.py]

Filter out certain TCRs based on a criterion from DeWitt_2018 (878 tcrs are filtered out by matching (v_family, amino_acid) tuple with those in a file from DeWitt 2018), and write out public tcr dict in csv format. There are 69,121,789 unique TCRs left and 8,739,207 unique public TCRs left. The output file is too big for github and is not included here.


## Compute TCR-HLA association p-values

[t6_reformat_HLA_v2_file.py]

reformat file of HLA existence among 666 HIP subjects to prepare for computing association p-values.

[t7_get_p_values.R]

compute p-value for association between public TCRs and HLAs using one-sided Fisher's exact test. This step took around 20 hours when using 36 cores on server.

## Get p-value cutoff under FDR 0.05

[t8_sort_pvalues.R]

sort true p-values.

[t9_get_p_values_permuted.R]

compute p-values under 20 permutations. The job was submitted using [t9_get_p_values_permuted.sh]. Each permutation took around 20 hours when using 36 cores on server.

[t10_get_counts.R]

under each permutation, for each true p-value, get the number of false positive under each true p-value cutoff. The binary search part in each run took around 50 minutes on server.

[t11_compute_pvalue_cutoff.R]

The new p-value cutoff is around 8.41189e-06 (more precisely, the p value on position 20,582 after all 1,878,929,505 are sorted).

## Get associated TCR (v allele-level) and HLA pairs

The file below runs interatively on Hutch cluster and writes out a csv file of the 20,582 associated pairs:

[t12_get_asso_pairs_v_allele.R]

Output file:

"../../data/intermediate_data/t12_HLA_associated_TCR_v_alleles.csv"

Based on the output file above, the following jupyter notebook

[t12_get_pairs_related_to_HLA_I_II.ipynb]

writes out five csv files:

The 6423 pairs related to HLA-I alleles as in HLA v2 of DeWitt_2018 (this version contains 85 HLA-I alleles, one in v1 is missing in v2)

"../../data/intermediate_data/t12_HLA_I_associated_TCR_v_alleles.csv"

The 14159 pairs related to HLA-II alleles as in HLA v2 of DeWitt_2018 (this version contains 125 simple HLA-II, DRB1, DPAB, or DQAB)

"../../data/intermediate_data/t12_HLA_II_associated_TCR_v_alleles.csv"

The 17281 pairs related to HLA-II alleles, similar as above, except that the 5 haplotypes are separated into 10 HLA-II

"../../data/intermediate_data/t12_HLA_II_expand_associated_TCR_v_alleles.csv"

The 11037 pairs for HLA-II alleles in the format of HLA v2, involving only the 125 simple ones, DRB1, DPAB, or DQAB, excluding the 5 haplotypes)

"../../data/intermediate_data/t12_HLA_II_no_haplotype_associated_TCR_v_alleles.csv"

The 6244 pairs involving the 5 haplotypes separated into 10 HLA-II

"../../data/intermediate_data/t12_HLA_II_5haplotypes_to_10pairs_associated_TCR_v_alleles.csv"

Among these, the first and the fourth file are used for later steps.


# Prepare Emerson data for building DePTH models


<br />  
