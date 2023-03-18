
# compute the sequence score matrix based on blosum62 using python
# load the results and compute distance based on the formula in 
# NetMHCpan paper

hla_filepath = "../data/Liu_2019/liu_2019_hla_ii_format.csv"
df_hla = read.csv(hla_filepath, header = TRUE)

df_blosum62_X_score = 
  read.csv("../results/st22_blosum62_X_seq_align_score_HLA_II.csv", header = TRUE)
            
# the column order of HLA matches that from HLA file
simi_mat_blosum62_X = matrix(NA, ncol = nrow(df_hla)[1], nrow = nrow(df_hla)[1])

for (i in 1:nrow(df_blosum62_X_score)){
  for (j in 1:nrow(df_blosum62_X_score)){
    simi_mat_blosum62_X[i, j] = 
      df_blosum62_X_score[i,j]/sqrt(df_blosum62_X_score[i,i]*df_blosum62_X_score[j,j])
  }
}

summary(c(simi_mat_blosum62_X))

dist_mat_blosum62_X = 1 - simi_mat_blosum62_X

summary(c(dist_mat_blosum62_X))

colnames(dist_mat_blosum62_X) = df_hla$hla

write.csv(dist_mat_blosum62_X, 
          file = "../results/st22_blosum62_HLA_II_dist.csv", 
          row.names = FALSE)

sessionInfo()
q(save="no")
