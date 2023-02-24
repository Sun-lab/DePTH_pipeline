# roc curves for HLA-I, 
# DePTH trained on Emerson 2017 and testted on Emerson 2017
# random forest trained on Emerson 2017 and tested on Emerson 2017

library(ggplot2)
library(pROC)


# load scores from DePTH

pos_depth_file = "HLA_I_all_match_test_pos_20.csv"
neg_depth_file = "HLA_I_all_match_test_neg_20.csv"

df_pos_depth = read.csv(paste0("../results/ensemble_scores/HLA_I_all_match/", pos_depth_file), header = TRUE)
df_neg_depth = read.csv(paste0("../results/ensemble_scores/HLA_I_all_match/", neg_depth_file), header = TRUE)

dim(df_pos_depth)
dim(df_neg_depth)

auc(c(rep(1, dim(df_pos_depth)[1]), rep(0, dim(df_neg_depth)[1])), c(df_pos_depth$ave, df_neg_depth$ave))


# load scores from random forest

pos_random_forest_file = "st3_test_pos_random_forest.csv"
neg_random_forest_file = "st3_test_neg_random_forest.csv"

df_pos_random_forest = read.csv(paste0("../results/", pos_random_forest_file), header = TRUE)
df_neg_random_forest = read.csv(paste0("../results/", neg_random_forest_file), header = TRUE)

dim(df_pos_random_forest)
dim(df_neg_random_forest)

auc(c(rep(1, dim(df_pos_random_forest)[1]), rep(0, dim(df_neg_random_forest)[1])), 
    c(df_pos_random_forest$random_forest_score, df_neg_random_forest$random_forest_score))


#define object to plot
rocobj_depth <- roc(c(rep(1, dim(df_pos_depth)[1]), rep(0, dim(df_neg_depth)[1])), 
              c(df_pos_depth$ave, df_neg_depth$ave))
rocobj_random_forest <- roc(c(rep(1, dim(df_pos_random_forest)[1]), rep(0, dim(df_neg_random_forest)[1])), 
                c(df_pos_random_forest$random_forest_score, df_neg_random_forest$random_forest_score))

rocs <- list()
rocs[[paste0("DePTH,              AUC 0.82")]] = rocobj_depth
rocs[[paste0("Random Forest, AUC 0.77")]] = rocobj_random_forest


p1 <- ggroc(rocs, aes=c("linetype", "color"), legacy.axes = TRUE) +
      geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1),
                   color="grey", linetype="dashed") + 
      theme(legend.position = c(0.5, 0.2)) + 
      theme(legend.title=element_blank()) + 
      theme_classic()

pdf(file = "../figures/depth_draft/step4_I_random_forest_roc_curve.pdf", width = 4.5, height = 2)
print(p1)
dev.off()



sessionInfo()
q(save="no")
