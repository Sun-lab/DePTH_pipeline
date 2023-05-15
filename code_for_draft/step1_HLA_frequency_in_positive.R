
library(ggplot2)

pos_file = "../data/intermediate_data/t12_HLA_I_associated_TCR_v_alleles.csv"
df_pos = read.csv(pos_file, header = TRUE)
dim(df_pos)
head(df_pos)

df_freq = data.frame(table(df_pos$hla_allele))
df_freq[order(df_freq$Freq, decreasing = TRUE), ]

p1 <- ggplot(df_freq, aes(x=Freq))+
      geom_histogram(color="darkblue", fill="lightblue") + 
      theme_classic() + xlab("freq among pos pairs")

pdf(file = "../figures/depth_draft/step1_freq_among_positive_pairs.pdf", width = 2.6, height = 2)
print(p1)
dev.off()