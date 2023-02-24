library(ggplot2)

df = read.csv("../results/st14_allelewise_auc_single.csv", header = TRUE)

df_reformat = data.frame(hla = c(df$hla, df$hla), 
                         auc = c(df$depth_mcpas, df$glazer_server), 
                         model = c(rep("DePTH McPAS (single)", 10), rep("CLAIRE", 10)))

df_reformat$model <- factor(df_reformat$model, levels = c("DePTH McPAS (single)", "CLAIRE"))

g2 = ggplot(data = df_reformat, aes(x = hla, y = auc, fill = model)) +
     geom_col(position = "dodge")  + theme_classic() + ylim(0, 1) + 
     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
     labs(x = "HLA", y = "AUC")

pdf("../figures/depth_draft/supp3_allelewise_auc_depth_mcpas_vs_claire_barchart_single.pdf", width = 4.8, height = 2.5)
print(g2)
dev.off()



sessionInfo()
q(save="no")