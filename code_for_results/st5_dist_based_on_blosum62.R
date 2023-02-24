
# compute the sequence score matrix based on blosum62 using python
# load the results and compute distance based on the formula in 
# NetMHCpan paper

HLA_filepath = "../data/HLA_I_pseudo_40.csv"
df_HLA = read.csv(HLA_filepath, header = TRUE)

df_blosum62_X_score = 
  read.csv("../results/st5_blosum62_X_seq_align_score_HLA_I.csv", header = TRUE)

# the column order of HLA matches that from HLA pseudo sequence file
simi_mat_blosum62_X = matrix(NA, ncol = nrow(df_HLA)[1], nrow = nrow(df_HLA)[1])

for (i in 1:nrow(df_blosum62_X_score)){
  for (j in 1:nrow(df_blosum62_X_score)){
    simi_mat_blosum62_X[i, j] = 
      df_blosum62_X_score[i,j]/sqrt(df_blosum62_X_score[i,i]*df_blosum62_X_score[j,j])
  }
}

summary(c(simi_mat_blosum62_X))

dist_mat_blosum62_X = 1 - simi_mat_blosum62_X

summary(c(dist_mat_blosum62_X))

colnames(dist_mat_blosum62_X) = df_HLA$hla

write.csv(dist_mat_blosum62_X, 
          file = "../results/st5_blosum62_HLA_I_dist.csv", 
          row.names = FALSE)

sessionInfo()
q(save="no")
