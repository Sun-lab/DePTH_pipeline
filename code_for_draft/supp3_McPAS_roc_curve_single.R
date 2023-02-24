
# roc curve for DePTH trained on Emerson 2017 and testted on Emerson 2017

library(ggplot2)
library(pROC)

# load scores from DePTH trained on McPAS training data and prediction on McPAS testing data
# (only with HLA-I alleles that we have pseudo sequence information for)

DePTH_neg_score_file = "HLA_I_full_mcpas_neg_test_5779_7821_6367/predicted_scores.csv"
DePTH_pos_score_file = "HLA_I_full_mcpas_pos_test_5779_7821_6367/predicted_scores.csv"

df_DePTH_neg_scores = read.csv(paste0("../results/predicted_scores/HLA_I_full_mcpas_ensemble/", DePTH_neg_score_file), header=TRUE)
df_DePTH_pos_scores = read.csv(paste0("../results/predicted_scores/HLA_I_full_mcpas_ensemble/", DePTH_pos_score_file), header=TRUE)

dim(df_DePTH_neg_scores)
dim(df_DePTH_pos_scores)

auc(c(rep(1, dim(df_DePTH_pos_scores)[1]), rep(0, dim(df_DePTH_neg_scores)[1])), 
    c(df_DePTH_pos_scores$score, df_DePTH_neg_scores$score))

# load scores from CLAIRE on server on McPAS test data
# (only with HLA-I alleles that we have pseudo sequence information for)

df_CLAIRE_pos_scores = 
  read.csv("../results/CLAIRE_mcpas/st2_Glazer_2022_online_test_kept_pos_scores.csv", header=TRUE)
df_CLAIRE_neg_scores = 
  read.csv("../results/CLAIRE_mcpas/st2_Glazer_2022_online_test_kept_neg_scores.csv", header=TRUE)

dim(df_CLAIRE_pos_scores)
dim(df_CLAIRE_neg_scores)

auc(c(rep(1, dim(df_CLAIRE_pos_scores)[1]), rep(0, dim(df_CLAIRE_neg_scores)[1])), 
    c(df_CLAIRE_pos_scores$score, df_CLAIRE_neg_scores$score))

#define object to plot
rocobj_i <- roc(c(rep(1, dim(df_DePTH_pos_scores)[1]), rep(0, dim(df_DePTH_neg_scores)[1])), 
              c(df_DePTH_pos_scores$score, df_DePTH_neg_scores$score))
rocobj_ii <- roc(c(rep(1, dim(df_CLAIRE_pos_scores)[1]), rep(0, dim(df_CLAIRE_neg_scores)[1])), 
                c(df_CLAIRE_pos_scores$score, df_CLAIRE_neg_scores$score))

rocs <- list()
rocs[[paste0("DePTH McPAS (single), AUC 0.76")]] = rocobj_i
rocs[[paste0("CLAIRE,                         AUC 0.78")]] = rocobj_ii


p1 <- ggroc(rocs, aes=c("linetype", "color"), legacy.axes = TRUE) +
      geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1),
                   color="grey", linetype="dashed") + 
      theme(legend.position = c(0.5, 0.2)) + 
      theme(legend.title=element_blank()) + 
      theme_classic()

pdf(file = "../figures/depth_draft/supp3_McPAS_roc_curve_single.pdf", width = 4.6, height = 2)
print(p1)
dev.off()

sessionInfo()
q(save="no")
