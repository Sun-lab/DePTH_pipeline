# cox regression for the association between survival outcome
# and HLA-I heterozygosity metrics
# based on Chowell 2018 subjects with HLA-I alleles that each either is or
# can be matched to one of the 85 Emerson HLA-I alleles with pseudo sequences

library(readxl)
library(survival)
library(survminer)
library(ggplot2)
library(ggpubr)
library(stringr)
theme_set(theme_bw())


df_2018_kept = read.csv("../results/st11_chowell_2018_depth_glazer_aa_set_breadth.csv", header = TRUE)
dim(df_2018_kept)
colSums(is.na(df_2018_kept))

table(df_2018_kept$gender)
table(df_2018_kept$drug_class)
table(df_2018_kept$cancer_type)

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

# change "" items in stage_m to NA
table(df_2018_kept$stage_m, useNA = "ifany")
na_stage_m = df_2018_kept$stage_m
na_stage_m[which(na_stage_m=="")] = NA
df_2018_kept$stage_m = na_stage_m
table(df_2018_kept$stage_m, useNA = "ifany")

# change "" items in stage to NA
table(df_2018_kept$stage, useNA = "ifany")
na_stage = df_2018_kept$stage
na_stage[which(na_stage=="")] = NA
df_2018_kept$stage = na_stage
table(df_2018_kept$stage, useNA = "ifany")

# neither stage_m nor stage is used because they both have NA for 
# a large number of subjects

# pad df_2018_kept mutcnt to prepare for log transformation

min(df_2018_kept$mutcnt[which(df_2018_kept$mutcnt>0)])

df_2018_kept$log_mutcnt = log10(df_2018_kept$mutcnt + 0.01)
summary(df_2018_kept$log_mutcnt)


# Explore association between different heterozygosity metrics and survival outcome

# survival analysis 


metrics_vec = c("homozygous", "depth_ave", "glazer_ave", "aa_ave", 
                "depth_set_ave", "glazer_set_ave", 
                "depth_breadth", "glazer_breadth")

rename_vec = c("homozygous", "mean_DePTH_cor", "mean_CLAIRE_cor", "mean_AA", 
               "mean_DePTH_set", "mean_CLAIRE_set", 
               "DePTH_breadth", "CLAIRE_breadth")
  
cur_matrix = matrix(NA, nrow = length(metrics_vec), 2)

# first, without covariates

col_i = 1

for (i in 1:length(metrics_vec)){
  x = metrics_vec[i]
  y = "Surv(os_months, os_event)"
  form = as.formula(paste(y, "~", x))
  cur_fit = coxph(form, data = df_2018_kept)
  cur_matrix[i, col_i] = summary(cur_fit)$coefficients[1, 5]
}

# next, with all covariates

col_i = 2

for (i in 1:length(metrics_vec)){
  
  x = metrics_vec[i]
  
  y = "Surv(os_months, os_event)"

  form = as.formula(paste(y, "~", x, "+ age_group + log_mutcnt + drug_class"))
  second_column = "age+log_mutcnt+drug_class"

  cur_fit = coxph(form, data = df_2018_kept)
  cur_matrix[i, col_i] = summary(cur_fit)$coefficients[1, 5]

}


df_cur_pvalues = as.data.frame(cur_matrix)
colnames(df_cur_pvalues) = c("no covariate", second_column)
rownames(df_cur_pvalues) = rename_vec

write.csv(df_cur_pvalues, 
          file = paste0("../results/st13_chowell_2018_survival_outcome_cox.csv"), 
          row.names = TRUE)  


  





gc()

sessionInfo()

q(save="no")


