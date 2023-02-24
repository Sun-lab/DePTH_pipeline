# this code does survival curves for clusters from step8_dist_cluster.ipynb
# with method = wald and 
# (1) height cutoff at 1.6 for DePTH set-based distance matrix



library(survival)
library(survminer)


# (1) height cutoffs 1.6 for DePTH set-based distance matrix

set_file = "../results/depth_draft/step8_chowell_2018_depth_glazer_aa_set_breadth_replace_depth_set_dist_hc_wald_1_1pt6.csv"
set_data = read.csv(set_file, header = TRUE)

set_data$cluster_cut1pt6 = as.character(set_data$cluster_cut1pt6)
table(set_data$cluster_cut1pt6)
table(set_data$reference, set_data$cluster_cut1pt6)



df_cut1pt6 = data.frame(group = set_data$cluster_cut1pt6, 
                        OS_Months = set_data$os_months, 
                        OS_Event = set_data$os_event)

cur_fit_cut1pt6 = coxph(Surv(OS_Months, OS_Event) ~ group, data = df_cut1pt6)
cur_fit_cut1pt6
# 0.039

fit_cut1pt6 <- survfit(Surv(OS_Months, OS_Event) ~ group, data=df_cut1pt6)

pdf("../figures/depth_draft/step9_chowell2018_depth_set_dist_hc_cut1pt6_survival_curves.pdf", width =4.2, height = 2.2)
print(ggsurvplot(fit_cut1pt6, data = df_cut1pt6,  
                 legend.labs = c("group 0", "group 1"), 
                 legend.title="", legend = "right"))
dev.off()

pdf("../figures/depth_draft/step9_chowell2018_depth_set_dist_hc_cut1pt6_survival_curves_n.pdf", width =4.8, height = 2.2)
print(ggsurvplot(fit_cut1pt6, data = df_cut1pt6,  
                 legend.labs = c("group 0 (n=909)", "group 1 (n=534)"), 
                 legend.title="", legend = "right"))
dev.off()




sessionInfo()

q(save="no")

