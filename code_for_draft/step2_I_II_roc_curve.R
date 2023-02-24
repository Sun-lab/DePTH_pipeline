# roc curve for DePTH trained on Emerson 2017 and testted on Emerson 2017

library(ggplot2)
library(pROC)


# load scores for hla_i

pos_i_file = "HLA_I_all_match_test_pos_20.csv"
neg_i_file = "HLA_I_all_match_test_neg_20.csv"

df_pos_i = read.csv(paste0("../results/ensemble_scores/HLA_I_all_match/", pos_i_file), header = TRUE)
df_neg_i = read.csv(paste0("../results/ensemble_scores/HLA_I_all_match/", neg_i_file), header = TRUE)

dim(df_pos_i)
dim(df_neg_i)

auc(c(rep(1, dim(df_pos_i)[1]), rep(0, dim(df_neg_i)[1])), c(df_pos_i$ave, df_neg_i$ave))


# load scores for hla_ii

pos_ii_file = "HLA_II_all_match_test_pos_20.csv"
neg_ii_file = "HLA_II_all_match_test_neg_20.csv"

df_pos_ii = read.csv(paste0("../results/ensemble_scores/HLA_II_all_match/", pos_ii_file), header = TRUE)
df_neg_ii = read.csv(paste0("../results/ensemble_scores/HLA_II_all_match/", neg_ii_file), header = TRUE)

dim(df_pos_ii)
dim(df_neg_ii)

auc(c(rep(1, dim(df_pos_ii)[1]), rep(0, dim(df_neg_ii)[1])), c(df_pos_ii$ave, df_neg_ii$ave))


#define object to plot
rocobj_i <- roc(c(rep(1, dim(df_pos_i)[1]), rep(0, dim(df_neg_i)[1])), 
              c(df_pos_i$ave, df_neg_i$ave))
rocobj_ii <- roc(c(rep(1, dim(df_pos_ii)[1]), rep(0, dim(df_neg_ii)[1])), 
                c(df_pos_ii$ave, df_neg_ii$ave))

rocs <- list()
rocs[[paste0("HLA", "\U002D", "I,  AUC 0.82")]] = rocobj_i
rocs[[paste0("HLA", "\U002D", "II, AUC 0.79")]] = rocobj_ii


p1 <- ggroc(rocs, aes=c("linetype", "color"), legacy.axes = TRUE) +
      geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1),
                   color="grey", linetype="dashed") + 
      theme(legend.position = c(0.5, 0.2)) + 
      theme(legend.title=element_blank()) + 
      theme_classic()

pdf(file = "../figures/depth_draft/step2_I_II_roc_curve.pdf", width = 4, height = 2)
print(p1)
dev.off()



sessionInfo()
q(save="no")
