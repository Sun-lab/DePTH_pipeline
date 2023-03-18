## Code files for getting results

st1:

  st1_settings_with_best_cv_valid_auc.ipynb

  combine the cross-validation average validation AUC from 60 hyperparameter settings into one file.

  st1_test_metrics.ipynb

  get the test AUC, recall and specificity from 20 models trained under 20 different sets of random seeds, under each data setting.


st2:

  st2_DePTH_all_match_full_mcpas.ipynb

  get the test AUC on McPAS based on scores given by DePTH trained on Emerson data, and the test AUC on Emerson data based on scores given by DePTH trained on McPAS data.

  st2_DePTH_full_mcpas_ensemble_AUCs.ipynb

  prepare result files needed for plotting the comparison between DePTH McPAS and CLAIRE (Glazer et al. 2022) model.


st3:

  st3_random_forest.ipynb

  run random forest on encoded Emerson data.


st4:

  st4_reshape_scores_on_zheng_2021_tcrs.ipynb

  convert predicted scores for pairs formed by 85 HLA-I alleles and 10008 Zheng 2021 TCRs to the format of matrix

  st4_dist_DePTH_zheng_2021_ensemble.R

  compute distance matrix for 85 HLA-I alleles based on correlation of predicted scores from DePTH on pairs with 10008 zheng 2021 TCRs


st5:

  st5_align_scores_based_on_blosum62.py

  st5_dist_based_on_blosum62.R

  compute distance matrix based on blosum_62


st6:

  st6_process_chowell_2018_data.ipynb

  process chowell 2018 data

  output all 146 HLA-I alleles that either belong to or can be matched to the 85 HLA-I alleles from Emerson data.

  output a file with 1443 subjects having HLA-I alleles belonging to the 146.


st7:

  st7_prepare_for_Glazer_2022_server_prediction.ipynb

  for the 146 HLA-I alleles from 1443 Chowell 2018 kept subjects, find possible matching to those involved in CLAIRE mcpas train+valid+test pairs and output a file for the correspondence.

  prepare files for getting prediction on pairs between 146 HLAs and 10008 zheng 2021 TCRs from CLAIRE(glazer 2022) server.


st8:

  st8_convert_Glazer_2022_server_prediction.ipynb

  convert prediction from CLAIRE (Glazer 2022) server on pairs between 146 HLAs and 10008 Zheng 2021 TCRs to the format of matrix


st9:

  st9_dist_Glazer_2022_server_matched_zheng_2021.R

  compute distance matrix among 146 HLA-I alleles from Chowell 2018 data based on predicted scores from CLAIRE (Glazer 2022) server
  on pairs with 10008 zheng 2021 TCRs.


st10:

  st10_chowell_2018_subject_level_dist_and_between_subject_ot_dist.ipynb

  for Chowell 2018 subjects, compute part of the HLA heterozygosity metrics, including mean_DePTH_cor, mean_CLAIRE_cor and mean_AA for subject-level metrics, and dist_DePTH_cor, dist_CLAIRE_cor and dist_AA for between-subject distances.


st11:

  st11_chowell_2018_tcr_set_distance.ipynb

  for Chowell 2018 subjects, compute part of the HLA heterozygosity metrics, including mean_DePTH_set, mean_CLAIRE_set, DePTH_breadth and CLAIRE_breadth for subject-level metrics, and dist_DePTH_breadth and dist_CLAIRE_breadth for between-subject distances.


st12:

   st12_chowell_2018_MiRKAT-S_test.R

   for Chowell 2018 subjects, run MiRKAT-S without or with covariates


st13:

   st13_chowell_2018_survival_outcome_cox.R

   for Chowell 2018 subjects, run cox regression without or with covariates


Flow chart for files related to MiRKAT-S and cox regression on chowell 2018 data:

```bash
                          st4 (DePTH) -->
                          st5 (AA)    -->| --> st10(chowell2018 1st part of metrics)  -> st12  MiRKAT-S
  st6 --> st7 --> st8 --> st9 (CLAIRE)-->      st11(chowell2018 2nd part of metrics)     st13  cox regression
   |---------------------------------------------------------------------------------------^
```

st14:

   st14_DePTH_full_mcpas_single_AUCs.ipynb

   prepare result files needed for plotting the comparison between DePTH McPAS (one single model, instead of ensemble) and CLAIRE (Glazer et al. 2022) model.

st15:

   st15_liu_2019_hla_i_identification.ipynb

   process Liu 2019 data

   extract 102 HLA-I alleles that either belong to or can be matched to the 85 HLA-I alleles from Emerson data.

   output a file with 143 subjects having HLA-I alleles belonging to 102 HLA-I alleles.

st16:

   st16_liu_2019_subject_level_dist.ipynb

   for Liu 2019 subjects, compute part of the HLA heterozygosity metrics, including mean_DePTH_cor and mean_AA for subject-level metrics.

st17:

   st17_liu_2019_tcr_set_distance.ipynb

   for Liu 2019 subjects, compute part of the HLA heterozygosity metrics, including mean_DePTH_set and DePTH_breadth for subject-level metrics.



Flow chart for files related to liu2019 HLA-I data processing

```bash
 st4 (DePTH) -->      
 st5 (AA)    -->| --> st16(liu2019 1st part of metrics)
 st15        -->      st17(liu2019 2nd part of metrics)
```

st18:

 st18_liu_2019_hla_ii_identification.ipynb

 process Liu 2019 data

 extract HLA-II components that either belong to or can be matched to the HLA-II alleles from Emerson data.

 output a file with 120 subjects having HLA-II components appearing in Emerson data.

st19:

 st19_liu_2019_extract_hla_ii_pairs.ipynb

 for kept Liu 2019 data, put DQA&DQB, DPA&DPB into combinations following the format of HLA-II alleles in Emerson data.

 output the list of unique HLA-II alleles.

st20:

 st20_reshape_scores_on_zheng_2021_ii_tcrs.ipynb

 reshape the scores got from DePTH for HLA-II alleles and Zheng 2021 CD4 positive TCRs.

st21:

 st21_dist_DePTH_zheng_2021_ii_ensemble.R

 compute distance matrix for kept Liu 2019 HLA-II alleles based on correlation of predicted scores from DePTH on pairs with zheng 2021 CD4 TCRs

st22:

 st22_align_hla_ii_scores_based_on_blosum62.py

 st22_hla_ii_dist_based_on_blosum62.R

 compute distance matrix for kept Liu 2019 HLA-II alleles based on blosum_62

st23:

 st23_liu_2019_hla_ii_subject_level_dist.ipynb

 for kept Liu 2019 subjects, compute part of the HLA heterozygosity metrics, including mean_DePTH_cor and mean_AA for subject-level metrics.

st24:

 st24_liu_2019_hla_ii_tcr_set_distance.ipynb

 for kept Liu 2019 subjects, compute part of the HLA heterozygosity metrics, including mean_DePTH_set and DePTH_breadth.



```bash  
                   |-------------------> st24
st18 --> st19 --> st20 --> st21 -->| --> st23
          |-------------> st22 --->|
```

st23 and st24 are followed by the code for plotting supplementary figures 4&5.
