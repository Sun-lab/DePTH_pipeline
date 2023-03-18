# this code computes distance matrices for HLAs


HLA_filepath = "../data/Liu_2019/liu_2019_hla_ii_format.csv"
df_HLA = read.csv(HLA_filepath, header = TRUE)
dim(df_HLA)

# distance based on positive TCRs from single cell data

pos_filename = "st20_HLA_II_zheng_2021_ensemble_reshape_20.csv" 
df_pos = read.csv(paste0("../results/", pos_filename), header = TRUE)
dim(df_pos)
head(df_pos)

stopifnot(dim(df_HLA)[1] == dim(df_pos)[2])

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
          file = "../results/st21_zheng_2021_pos_HLA_II_dist_ensemble_20.csv", 
          row.names = FALSE)


gc()

sessionInfo()
q(save="no")