
library(ggplot2)

pos_file = "../data/step77_HLA_I_associated_TCR_v_alleles.csv"
df_pos = read.csv(pos_file, header = TRUE)
dim(df_pos)
head(df_pos)

df_freq = data.frame(table(df_pos$hla_allele))

p1 <- ggplot(df_freq, aes(x=Freq))+
      geom_histogram(color="darkblue", fill="lightblue") + 
      theme_classic() + xlab("freq among pos pairs")

pdf(file = "../figures/depth_draft/step2_freq_among_positive_pairs.pdf", width = 2.6, height = 2)
print(p1)
dev.off()