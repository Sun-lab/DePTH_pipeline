
# Pipeline for DePTH

Folder [process_emerson_data]

Code files for processing Emerson data, find pseudo sequences for HLAs and getting p-value cutoff for association.

Folder [prepare_emerson_data_sets_for_DePTH]

Construct positive and negative TCR-HLA pairs for building DePTH models based on processed Emerson data.

**Szeto 2020:**

[a1_prepare_szeto_2020_pairs.ipynb]

Process Szeto 2020 solved structures to get the 54 pairs for DePTH and Random Forest to make prediction on.

**Glazer 2020 McPAS:**

[a2_prepare_glazer_2022_full_mcpas.ipynb]

Process McPAS data used by Glazer et al. 2022 for training CLAIRE model, to generate pairs based on McPAS data to build DePTH McPAS models.

**Zheng 2021:**

[a3_process_CD8T.Rmd]

[a4_extract_TCR_CD8_for_DePTH.Rmd]

[a5_convert_zheng_2021_CD8_format.ipynb]

[a6_prepare_zheng_2021_CD8_pairs.ipynb]

Extract potentially cancer-related CD8 TCRs based on Zheng 2021 data and gene signatures from other studies, and prepare the input data file based on 85 HLA-I alleles from Emerson data and 10008 processed potentially cancer-related CD8 TCRs for DePTH to make prediction on.

[a3_process_CD4T.Rmd]

[a4_extract_TCR_CD4_for_DePTH.Rmd]

[a7_convert_zheng_2021_CD4_format.ipynb]

[a8_prepare_zheng_2021_CD4_liu_2019_pairs.ipynb]

Extract potentially cancer-related CD4 TCRs based on Zheng 2021 data and gene signatures from other studies, and prepare the input data file based on 141 HLA-II alleles from the intersection of processed Liu 2019 data and Emerson data and 6547 processed potentially cancer-related CD4 TCRs for DePTH to make prediction on.


<br />  
