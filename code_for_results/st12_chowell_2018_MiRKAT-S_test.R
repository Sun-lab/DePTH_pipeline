

# test association with survival outcome, on Chowell 2018 subjects
#  -- DePTH on Zheng 2021 pos tcrs, cor (dist_DePTH_cor)
#  -- Glazer 2022 on Zheng 2021 pos tcrs, cor (dist_CLAIRE_cor)
#  -- aa blosum62 (dist_AA)
#  -- DePTH on Zheng 2021, set (dist_DePTH_breadth)
#  -- Glazer 2022 on Zheng 2021, set (dist_CLAIRE_breadth)

library(MiRKAT)

# load distance matrices

df_depth_cor = read.csv("../results/st10_chowell_2018_kept_depth_ot_mat.csv", header = TRUE)
df_glazer_cor = read.csv("../results/st10_chowell_2018_kept_glazer_ot_mat.csv", header = TRUE)
df_aa = read.csv("../results/st10_chowell_2018_kept_blosum62_ot_mat.csv", header = TRUE)
df_depth_set = read.csv("../results/st11_chowell_2018_kept_depth_set_mat.csv", header = TRUE)
df_glazer_set = read.csv("../results/st11_chowell_2018_kept_glazer_set_mat.csv", header = TRUE)

dim(df_depth_cor)
dim(df_glazer_cor)
dim(df_aa)
dim(df_depth_set)
dim(df_glazer_set)

# check whether the matrices are symmetric
mat_depth_cor = as.matrix(df_depth_cor)
mat_depth_cor_diff = abs(mat_depth_cor - t(mat_depth_cor))
max(mat_depth_cor_diff)

mat_glazer_cor = as.matrix(df_glazer_cor)
mat_glazer_cor_diff = abs(mat_glazer_cor - t(mat_glazer_cor))
max(mat_glazer_cor_diff)

mat_aa = as.matrix(df_aa)
mat_aa_diff = abs(mat_aa - t(mat_aa))
max(mat_aa_diff)

mat_depth_set = as.matrix(df_depth_set)
mat_depth_set_diff = abs(mat_depth_set - t(mat_depth_set))
max(mat_depth_set_diff)

mat_glazer_set = as.matrix(df_glazer_set)
mat_glazer_set_diff = abs(mat_glazer_set - t(mat_glazer_set))
max(mat_glazer_set_diff)



df_2018_kept = read.csv("../results/st11_chowell_2018_depth_glazer_aa_set_breadth.csv", header = TRUE)
dim(df_2018_kept)
colSums(is.na(df_2018_kept))

table(df_2018_kept$gender, useNA="ifany")
table(df_2018_kept$drug_class, useNA="ifany")
table(df_2018_kept$cancer_type, useNA="ifany")

reset_drug_class = df_2018_kept$drug_class
reset_drug_class[which(reset_drug_class%in%c("PD-1/PDL-1", "PD-L1"))] = "PD-1"
reset_drug_class[which(reset_drug_class%in%c("CTLA-4 + PD-1", "PD-L1 + CTLA-4"))] = "Combo"

df_2018_kept$drug_class = reset_drug_class
table(df_2018_kept$drug_class)

table(df_2018_kept$reference)

df_2018_kept$age_group[which(df_2018_kept$age_group == "<30")] = 1
df_2018_kept$age_group[which(df_2018_kept$age_group == "31-50")] = 2
df_2018_kept$age_group[which(df_2018_kept$age_group == "50-60")] = 3
df_2018_kept$age_group[which(df_2018_kept$age_group == "61-70")] = 4
df_2018_kept$age_group[which(df_2018_kept$age_group == ">71")] = 5

df_2018_kept$age_group = as.numeric(df_2018_kept$age_group)


table(df_2018_kept$stage_m, useNA = "ifany")
table(df_2018_kept$stage, useNA = "ifany")


# pad df_2018_kept mutcnt to prepare for log transformation

min(df_2018_kept$mutcnt[which(df_2018_kept$mutcnt>0)])

df_2018_kept$log_mutcnt = log10(df_2018_kept$mutcnt + 0.01)
summary(df_2018_kept$log_mutcnt)

colnames(df_2018_kept)


table(is.na(df_2018_kept$log_mutcnt), df_2018_kept$drug_class, useNA="ifany")

table(is.na(df_2018_kept$log_mutcnt), df_2018_kept$cancer_type, useNA="ifany")

table(df_2018_kept$stage_m, df_2018_kept$drug_class, useNA="ifany")


# age_group + log_mutcnt + drug_class
# because stage_m is missing a large number of subjects


# define a function to build feature matrix


process_mat_feature <- function(cur_mat, cur_df, cur_features){
    
    stopifnot(nrow(cur_mat)==nrow(cur_df))

    temp_df = cur_df[, cur_features]

    cur_rows = which(rowSums(is.na(temp_df))==0)

    output_mat = cur_mat[cur_rows, cur_rows]

    nonan_df = temp_df[cur_rows, ]
    output_df = cur_df[cur_rows, ]
    
    len_features = length(cur_features)

    output_features = c()


    for (i in 1:len_features){

        ft = cur_features[i]
        
        if (class(nonan_df[[ft]])=="character"){
          
          uni_values = unique(nonan_df[[ft]])

          if (length(uni_values)>1){

            for (j in 2:length(uni_values)){

                output_features = cbind(output_features, 
                                        as.numeric(nonan_df[[ft]]==uni_values[j]))
            }
          }
        }else{
          
          output_features = cbind(output_features, 
                                  nonan_df[[ft]])
        }
    }

   output_list = list()
   output_list[["mat"]] = output_mat
   output_list[["X"]] = output_features
   output_list[["surv"]] = output_df


   return(output_list)
    
}



mat_list = list()
mat_list[["dist_DePTH_cor"]] = as.matrix(df_depth_cor)
mat_list[["dist_CLAIRE_cor"]] = as.matrix(df_glazer_cor)
mat_list[["dist_AA"]] = as.matrix(df_aa)
mat_list[["dist_DePTH_breadth"]] = as.matrix(df_depth_set)
mat_list[["dist_CLAIRE_breadth"]] = as.matrix(df_glazer_set)

matrix_names = names(mat_list)

wo_vec = rep(NA, 5)
w_vec = rep(NA, 5)

feature_vec = c("age_group", "log_mutcnt", "drug_class")


for (i in 1:5){
  
  name_i = matrix_names[i]
  mat_i = mat_list[[name_i]]
  k_i = D2K(mat_i)
  
  wo_vec[i] = MiRKATS(Ks = k_i, 
                      obstime = df_2018_kept$os_months, 
                      delta = df_2018_kept$os_event)$p_values
  
}

wo_vec


for (i in 1:5){
    
  name_i = matrix_names[i]
  mat_i = mat_list[[name_i]]
  list_i = process_mat_feature(mat_i, df_2018_kept, feature_vec)
  k_i = D2K(list_i[["mat"]])

  w_vec[i] = MiRKATS(Ks = k_i, 
                     obstime = list_i[["surv"]]$os_months, 
                     delta = list_i[["surv"]]$os_event, 
                     X = list_i[["X"]])$p_values

}

w_vec



df_pvalues = data.frame(wo_covariate = wo_vec,
                        w_covariates = w_vec)

colnames(df_pvalues) = c("no_covariate", "age+log_mutcnt+drug_class")
row.names(df_pvalues) = matrix_names
  
df_pvalues

write.csv(df_pvalues, 
          file = "../results/st12_chowell_2018_MiRKAT-S_pvalues.csv", 
          row.names = TRUE)

gc()

sessionInfo()
q(save="no")