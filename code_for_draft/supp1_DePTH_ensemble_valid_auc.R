# from 100 models trained using different sets of random seeds
# each time, 
#     do permutation
#     for n = 1, ..., 100, get the average validation AUC from the first n models

  

library(ggplot2)
theme_set(theme_classic())
library(ggpubr)
library(pROC)


df_seeds = read.table("ensemble_seeds.txt", sep = "\t", header=FALSE)
dim(df_seeds)
head(df_seeds)

score_dir = "HLA_I_model_ensemble_prediction/"

valid_pos_list = c()
valid_neg_list = c()


for (i in 1:100){
  
  seed_1 = as.character(df_seeds$V1[i])
  seed_2 = as.character(df_seeds$V2[i])
  seed_3 = as.character(df_seeds$V3[i])
  
  seeds_string = paste0("_", seed_1, "_", seed_2, "_", seed_3)
  
  df_valid_pos = read.csv(paste0(score_dir, "HLA_I_pos_valid_load", seeds_string, 
                                 "/predicted_scores.csv"), 
                          header=TRUE)
  df_valid_neg = read.csv(paste0(score_dir, "HLA_I_neg_valid_load", seeds_string, 
                                 "/predicted_scores.csv"), 
                          header=TRUE)  

  
  valid_pos_list[[i]] = df_valid_pos$score
  valid_neg_list[[i]] = df_valid_neg$score
  
}

length(valid_pos_list)
length(valid_neg_list)

valid_pos_mat = matrix(unlist(valid_pos_list), ncol=100, byrow=FALSE)
valid_neg_mat = matrix(unlist(valid_neg_list), ncol=100, byrow=FALSE)

dim(valid_pos_mat)
dim(valid_neg_mat)

stopifnot(all(valid_pos_mat[1:5, 20] == valid_pos_list[[20]][1:5]))
stopifnot(all(valid_neg_mat[1:5, 30] == valid_neg_list[[30]][1:5]))



# Validation AUC
# AUC based on accumulative average, four shuffles
# In each round, first, shuffle all 100 models. 
# Then, record the validation AUC based on the average scores of
#   the first 1, 2, ..., 100 models. 

p_valid_auc_shuffled = list()

valid_labels = c(rep(1, nrow(valid_pos_mat)), rep(0, nrow(valid_neg_mat)))
  
set.seed(1629)

for (h in 1:4){
  
  shuffle_indexes = sample(100, 100, replace=FALSE)
  
  accum_valid_aucs_shuffled = rep(NA, 100)
  accum_valid_aucs_shuffled[1] = auc(valid_labels, 
                                     c(valid_pos_mat[, shuffle_indexes[1]], 
                                       valid_neg_mat[, shuffle_indexes[1]]))
  
  for (k in 2:100){
    
    cur_V1 = rowSums(valid_pos_mat[, shuffle_indexes[1:k]])/k
    cur_V2 = rowSums(valid_neg_mat[, shuffle_indexes[1:k]])/k
    
    accum_valid_aucs_shuffled[k] = auc(valid_labels, c(cur_V1, cur_V2))
  
  }

  df_accum_valid_shuffled = data.frame(n_model=1:100, 
                                       valid_auc=accum_valid_aucs_shuffled)
  
  cur_p = ggplot(df_accum_valid_shuffled, aes(x=n_model, y=valid_auc)) + geom_point() + 
          xlab("number of models averaged") + 
          ylab("validation AUC") + 
          ggtitle(paste0("Shuffle #", as.character(h)))
  
  p_valid_auc_shuffled[[h]] = cur_p

}

pdf(file = "supp1_validation_auc_ave_first_n_models.pdf", width = 6.6, height = 5.2)
ggarrange(plotlist=p_valid_auc_shuffled, ncol=2, nrow = 2)
dev.off()



sessionInfo()
q(save="no")
