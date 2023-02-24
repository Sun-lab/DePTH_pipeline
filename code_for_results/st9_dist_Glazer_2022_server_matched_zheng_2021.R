# this code computes distance matrices for HLAs, based on predicted scores
# from Glazer 2022 model on server on the pairs between Zheng 2021 tcrs and
# HLAs


# distance based on positive TCRs from single cell data

pos_filename = "st8_Glazer_2022_server_on_zheng_2021_pos_extended_hlas_146.csv" 
df_pos = read.csv(paste0("../results/", pos_filename), header = TRUE)
dim(df_pos)
head(df_pos)


corr_pos = matrix(NA, ncol = dim(df_pos)[2], nrow = dim(df_pos)[2])
corr_pos[dim(df_pos)[2], dim(df_pos)[2]] = 1

for (i in 1:(dim(df_pos)[2]-1)){
    corr_pos[i, i] = 1
  for (j in (i+1):(dim(df_pos)[2])){
    x = as.numeric(df_pos[, i])
    y = as.numeric(df_pos[, j])
    corr_pos[i, j] = cor(x, y, method = c("spearman"))
    corr_pos[j, i] = corr_pos[i, j]
  }
}

dist_pos = 1 - corr_pos

summary(c(dist_pos))

colnames(dist_pos) = colnames(df_pos)


write.csv(as.data.frame(dist_pos), 
          file = "../results/st9_Glazer_2022_server_on_Zheng_2021_pos_HLA_I_dist.csv", 
          row.names = FALSE)


gc()

sessionInfo()
q(save="no")