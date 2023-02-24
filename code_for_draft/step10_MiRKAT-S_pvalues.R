# this code draws barplot for the p-values from MiRKAT-S on Chowell 2018
# subjects after filtering out subjects with HLA-Is that cannot be matched 
# to any of the 85 ones from Emerson

library(ggplot2)
library(ggpubr)


result_dir = ("../results")
MiRKAT_S_2018_table = "st12_chowell_2018_MiRKAT-S_pvalues.csv"

df_mirkats_2018 = read.csv(paste0(result_dir, "/", MiRKAT_S_2018_table), header = TRUE)

colnames(df_mirkats_2018) = c("metric", "no_covariate", "w_covariate")

df_mirkats_2018$metric <- factor(df_mirkats_2018$metric, 
                                 levels = c("dist_AA", "dist_CLAIRE_cor", "dist_DePTH_cor", 
                                            "dist_CLAIRE_breadth", "dist_DePTH_breadth"))


p1 = ggplot(df_mirkats_2018, aes(x=metric, y=-log10(no_covariate), fill=metric)) +
     geom_bar(stat="identity")+theme_classic() + ylab("-log10(p value)") + 
     geom_hline(yintercept=-log10(0.05), color = "grey") + 
     theme(axis.text.x=element_blank(),
           axis.ticks.x=element_blank()) + 
     theme(legend.title=element_blank()) + 
     xlab("between-subject distance") + 
     ylim(0, 2.2) + ggtitle("Without covariates")



p2 = ggplot(df_mirkats_2018, aes(x=metric, y=-log10(w_covariate), fill=metric)) +
     geom_bar(stat="identity")+theme_classic() + ylab("-log10(p value)") + 
     geom_hline(yintercept=-log10(0.05), color = "grey") + 
     theme(axis.text.x=element_blank(),
           axis.ticks.x=element_blank()) + 
     theme(legend.title=element_blank()) + 
     xlab("between-subject distance") + 
     ylim(0, 2.2) + ggtitle("With covariates")


pdf(file = "../figures/depth_draft/step10_mirkats_chowell_2018.pdf", width = 6.6, height = 2)
ggarrange(p1, p2, nrow = 1, ncol = 2, 
          common.legend = TRUE, legend="right")
dev.off()


sessionInfo()
q(save="no")