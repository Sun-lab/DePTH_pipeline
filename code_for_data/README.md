
# Pipeline for DePTH

Folder [process_emerson_data]

Code files for processing Emerson data, find pseudo sequences for HLAs and getting p-value cutoff for association.

Folder [prepare_emerson_data_sets_for_DePTH]

Construct positive and negative TCR-HLA pairs for building DePTH models based on processed Emerson data.

[a1_prepare_szeto_2020_pairs.ipynb]

Process Szeto 2020 solved structures to get the 54 pairs for DePTH and Random Forest to make prediction on.

[a2_prepare_zheng_2021_pairs.ipynb]

Prepare the input data file based on 85 HLA-I alleles from Emerson data and 10008 processed potentially cancer-related Zheng 2021 TCRs for DePTH to make prediction on.

[a3_prepare_glazer_2022_full_mcpas.ipynb]

Process McPAS data used by Glazer et al. 2022 for training CLAIRE model, to generate pairs based on McPAS data to build DePTH McPAS models.




<br />  
