# this code gets the side-by-side violin plots based on the scores on 
# test positive pairs, test negative pairs, Szeto 2021 43 pairs

library(ggplot2)
library(ggpubr)
library(pROC)

theme_set(theme_classic())

results_dir = "../results/"

pos_random_forest_file = "st3_test_pos_random_forest.csv"
neg_random_forest_file = "st3_test_neg_random_forest.csv"

df_test_pos_sc = read.csv(paste0(results_dir, pos_random_forest_file), header = TRUE)
df_test_neg_sc = read.csv(paste0(results_dir, neg_random_forest_file), header = TRUE)

dim(df_test_pos_sc)
dim(df_test_pos_sc)


df_szeto_test_pos_sc = read.csv(paste0(results_dir, "st3_szeto_random_forest.csv"), header = TRUE)

dim(df_szeto_test_pos_sc)


auc(c(rep(1, nrow(df_szeto_test_pos_sc)), rep(0, nrow(df_test_neg_sc))), 
    c(df_szeto_test_pos_sc$random_forest_score, df_test_neg_sc$random_forest_score))

score_vec = c(df_test_pos_sc$random_forest_score, 
              df_test_neg_sc$random_forest_score, 
              df_szeto_test_pos_sc$random_forest_score)
data_set_vec = c(rep("test_pos", nrow(df_test_pos_sc)), 
                 rep("test_neg", nrow(df_test_neg_sc)), 
                 rep("Szeto_2020", nrow(df_szeto_test_pos_sc)))

df_scores = data.frame(dataset = data_set_vec, 
                       score = score_vec)


df_scores$dataset = factor(df_scores$dataset, levels=c("test_pos", "test_neg", "Szeto_2020"))

pdf(file = "../figures/depth_draft/step5_random_forest_violin_emerson_szeto.pdf", width = 3.5, height = 2)
ggplot(df_scores, aes(x=dataset, y=score, color=dataset)) + 
  xlab("data set") +
  ylab("predicted score") + 
  geom_violin(trim=FALSE) + 
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), 
        legend.title=element_blank())
#+ theme(legend.position="none")
dev.off()




sessionInfo()
q(save="no")
