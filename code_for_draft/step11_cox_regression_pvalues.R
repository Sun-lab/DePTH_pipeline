# this code draws barplot for the p-values from cox regression on 
# Chowell 2018 patients with HLA-Is that can be matched to one of
# the 85 HLA-I alleles from Emerson data

library(ggplot2)
library(ggpubr)

result_dir = ("../results")
cox_table = "st13_chowell_2018_survival_outcome_cox.csv"

df = read.csv(paste0(result_dir, "/", cox_table), header = TRUE)
df
colnames(df)[1] = "metric"

df$metric <- factor(df$metric, 
                   levels = c("homozygous", "mean_AA", 
                              "mean_CLAIRE_cor", "mean_DePTH_cor",
                              "mean_CLAIRE_set", "mean_DePTH_set", 
                              "CLAIRE_breadth", "DePTH_breadth"))


p1 = ggplot(df, aes(x=metric, y=-log10(no.covariate), fill=metric)) +
     geom_bar(stat="identity")+theme_classic() + ylab("-log10(p value)") + 
     geom_hline(yintercept=-log10(0.05), color = "grey") + 
     ggtitle("Without covariates") + 
     ylim(0, 2.35) + 
     theme(axis.text.x=element_blank(),
           axis.ticks.x=element_blank()) + 
     theme(legend.title=element_blank()) + 
     xlab("subject-level metric") 

p2 = ggplot(df, aes(x=metric, y=-log10(age.log_mutcnt.drug_class), fill=metric)) +
  geom_bar(stat="identity")+theme_classic() + ylab("-log10(p value)") + 
  geom_hline(yintercept=-log10(0.05), color = "grey") + 
  ggtitle("With covariates") + 
  ylim(0, 2.35) + 
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + 
  theme(legend.title=element_blank()) + 
  xlab("subject-level metric") 


pdf(file = "../figures/depth_draft/step11_cox_regression_chowell_2018.pdf", width = 7.6, height = 2.4)
ggarrange(p1, p2, nrow = 1, ncol = 2, 
          common.legend = TRUE, legend="right")
dev.off()







# p11 = ggplot(df_cox_2019_all, aes(x=metric, y=-log10(no.covariate), fill=metric)) +
#   geom_bar(stat="identity")+theme_classic() + ylab("-log10(p value)") + 
#   geom_hline(yintercept=-log10(0.05), color = "grey") + 
#   ggtitle("All") + 
#   ylim(0, 2.35) + 
#   theme(axis.text.x=element_blank(),
#         axis.ticks.x=element_blank())
# 
# pdf(file = "../figures/depth_draft/step17_figures_legend_right_main/step17_main_cox_chowell_2019_all.pdf", width = 5, height = 3)
# print(p11)
# dev.off()
# 
# p12 = ggplot(df_cox_2019_exome, aes(x=metric, y=-log10(no.covariate), fill=metric)) +
#   geom_bar(stat="identity")+theme_classic() + ylab("-log10(p value)") + 
#   geom_hline(yintercept=-log10(0.05), color = "grey") + 
#   ggtitle("Exome") + 
#   ylim(0, 2.35) + 
#   theme(axis.text.x=element_blank(),
#         axis.ticks.x=element_blank())
# 
# pdf(file = "../figures/depth_draft/step17_figures_legend_right_main/step17_main_cox_chowell_2019_exome.pdf", width = 5, height = 3)
# print(p12)
# dev.off()
# 
# p13 = ggplot(df_cox_2019_panel, aes(x=metric, y=-log10(no.covariate), fill=metric)) +
#   geom_bar(stat="identity")+theme_classic() + ylab("-log10(p value)") + 
#   geom_hline(yintercept=-log10(0.05), color = "grey") + 
#   ggtitle("Panel") + 
#   ylim(0, 2.35) + 
#   theme(axis.text.x=element_blank(),
#         axis.ticks.x=element_blank())
# 
# pdf(file = "../figures/depth_draft/step17_figures_legend_right_main/step17_main_cox_chowell_2019_panel.pdf", width = 5, height = 3)
# print(p13)
# dev.off()





sessionInfo()
q(save="no")