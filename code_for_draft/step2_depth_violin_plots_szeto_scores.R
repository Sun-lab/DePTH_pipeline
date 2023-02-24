# this code gets the side-by-side violin plots based on the scores on 
# test positive pairs, test negative pairs, Szeto 2021 54 pairs

library(ggplot2)
library(ggpubr)
library(pROC)

theme_set(theme_classic())

results_dir = "../results/ensemble_scores/"

res_emerson_szeto_folder = "HLA_I_all_match_szeto/"
res_emerson_emerson_folder = "HLA_I_all_match/"


df_emerson_emerson_test_pos_sc = read.csv(paste0(results_dir, res_emerson_emerson_folder,  
                                                 "HLA_I_all_match_test_pos_20.csv"), header = TRUE)

df_emerson_emerson_test_neg_sc = read.csv(paste0(results_dir, res_emerson_emerson_folder, 
                                                 "HLA_I_all_match_test_neg_20.csv"), header = TRUE)

dim(df_emerson_emerson_test_pos_sc)
dim(df_emerson_emerson_test_neg_sc)


df_emerson_szeto_test_pos_sc = read.csv(paste0(results_dir, res_emerson_szeto_folder, 
                                              "HLA_I_all_match_szeto_20.csv"), header = TRUE)

dim(df_emerson_szeto_test_pos_sc)


auc(c(rep(1, nrow(df_emerson_szeto_test_pos_sc)), rep(0, nrow(df_emerson_emerson_test_neg_sc))), 
    c(df_emerson_szeto_test_pos_sc$ave, df_emerson_emerson_test_neg_sc$ave))

emerson_score_vec = c(df_emerson_emerson_test_pos_sc$ave, 
                      df_emerson_emerson_test_neg_sc$ave, 
                      df_emerson_szeto_test_pos_sc$ave)
emerson_data_set_vec = c(rep("test_pos", nrow(df_emerson_emerson_test_pos_sc)), 
                         rep("test_neg", nrow(df_emerson_emerson_test_neg_sc)), 
                         rep("Szeto_2020", nrow(df_emerson_szeto_test_pos_sc)))

df_emerson_scores = data.frame(dataset = emerson_data_set_vec, 
                               score = emerson_score_vec)


df_emerson_scores$dataset = factor(df_emerson_scores$dataset, levels=c("test_pos", "test_neg", "Szeto_2020"))

pdf(file = "../figures/depth_draft/step2_depth_violin_emerson_szeto.pdf", width = 3.5, height = 2)
ggplot(df_emerson_scores, aes(x=dataset, y=score, color=dataset)) + 
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
