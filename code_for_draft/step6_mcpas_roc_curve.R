# roc curve for DePTH trained on mcpas and tested on mcpas
# in comparison with CLAIRE on mcpas

library(ggplot2)
library(pROC)


# load scores from DePTH trained on McPAS training data and prediction on McPAS testing data
# (only with HLA-I alleles that we have pseudo sequence information for)

DePTH_neg_score_file = "HLA_I_full_mcpas_test_neg_20.csv"
DePTH_pos_score_file = "HLA_I_full_mcpas_test_pos_20.csv"

df_DePTH_neg_scores = read.csv(paste0("../results/ensemble_scores/HLA_I_full_mcpas/", DePTH_neg_score_file), header=TRUE)
df_DePTH_pos_scores = read.csv(paste0("../results/ensemble_scores/HLA_I_full_mcpas/", DePTH_pos_score_file), header=TRUE)

auc(c(rep(1, dim(df_DePTH_pos_scores)[1]), rep(0, dim(df_DePTH_neg_scores)[1])), 
    c(df_DePTH_pos_scores$ave, df_DePTH_neg_scores$ave))


# load scores from CLAIRE on server on McPAS testing data
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
              c(df_DePTH_pos_scores$ave, df_DePTH_neg_scores$ave))
rocobj_ii <- roc(c(rep(1, dim(df_CLAIRE_pos_scores)[1]), rep(0, dim(df_CLAIRE_neg_scores)[1])), 
                c(df_CLAIRE_pos_scores$score, df_CLAIRE_neg_scores$score))

rocs <- list()
rocs[[paste0("DePTH McPAS, AUC 0.79")]] = rocobj_i
rocs[[paste0("CLAIRE,            AUC 0.78")]] = rocobj_ii


p1 <- ggroc(rocs, aes=c("linetype", "color"), legacy.axes = TRUE) +
      geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1),
                   color="grey", linetype="dashed") + 
      theme(legend.position = c(0.5, 0.2)) + 
      theme(legend.title=element_blank()) + 
      theme_classic()

pdf(file = "../figures/depth_draft/step6_mcpas_roc_curve.pdf", width = 4.45, height = 2)
print(p1)
dev.off()





sessionInfo()
q(save="no")
