# for HLA-B08:01, based on the patients with known existing or not information of this HLA
# plot histograms of frequency of positive TCRs and that of negative TCRs
# overlay the two histograms

library(ggplot2)

pos_file = "../data/intermediate_data/t12_HLA_I_associated_TCR_v_alleles.csv"
df_pos = read.csv(pos_file, header = TRUE)
dim(df_pos)
head(df_pos)

df_pos_B0801 = df_pos[which(df_pos$hla_allele=="HLA-B*08:01"), ]

freq_file = "../data/intermediate_data/HLA_I_all_match_prepare/tcr_freq_HLA_I_1.csv"
df_freq = read.csv(freq_file, header = TRUE)
dim(df_freq)
head(df_freq)



df_freq$full_tcr = paste0(df_freq$v_allele, ",", df_freq$amino_acids)

table(df_pos_B0801$tcr %in% df_freq$full_tcr)

col_group = rep(NA, dim(df_freq)[1])
col_group[which(df_freq$full_tcr%in%df_pos_B0801$tcr)] = "pos"
col_group[which(!(df_freq$full_tcr%in%df_pos_B0801$tcr))] = "neg"

df_freq$TCR_group = col_group
df_freq$TCR_group = factor(df_freq$TCR_group, levels = c("pos", "neg"))

head(df_freq)

p1 <- ggplot(df_freq,aes(x=log10(freq),fill=TCR_group))+
      geom_histogram(aes(y=..density..), 
                     alpha=0.5,position='identity',binwidth=0.2) + 
      theme_classic() + xlab("log10(population frequency)") +
      guides(fill=guide_legend(title="TCR group"))

pdf(file = "../figures/depth_draft/step1_tcr_population_frequency.pdf", width = 3.6, height = 2)
print(p1)
dev.off()

sessionInfo()
q(save="no")
