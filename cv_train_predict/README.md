
# Scripts and code files for cv, training and prediction

There are six situations of building DePTH models based on different data sets, HLA_I_all_match, HLA_II_all_match, HLA_I_leave1out_match (three situations included) and HLA_I_full_mcpas.

In each situation, first run cv script through 60 different hyperparameter settings. After combining validation AUCs from 60 settings into one file (using code as in code_for_results/st1_settings_with_best_cv_valid_auc.ipynb), the best hyperparameter setting is used to submit training script to train models using 20 different sets of random seeds. Eventually, for each test TCR-HLA pair, the average of predicted scores from the 20 trained models is taken as the final predicted score.

In the situation of HLA_I_all_match, 100 models using 100 different sets of random seeds are trained for the purpose of doing shuffling and choosing appropriate ensemble size (code file code_for_draft/supp1_DePTH_ensemble_valid_auc.R). In addition, predicted scores based on the scores of 20 trained model on processed pairs from Szeto 2020, pairs between 85 HLA-I alleles from Emerson data and 10008 TCRs from Zheng 2021 potentially cancer-related TCRs are also obtained. 






<br />  
