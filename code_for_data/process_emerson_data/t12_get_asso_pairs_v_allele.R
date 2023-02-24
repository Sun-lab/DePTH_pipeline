# this file writes out the associated (TCR, HLA) pairs based on p values 
# from v allele level TCR into a csv file of 
# the associated pairs and corresponding p values

data_dir = "../../data/intermediate_data/"

pvalues_big_list = readRDS(paste0(data_dir, "pvalues_big_list.rds"))
pvalues = unlist(pvalues_big_list)


# pvalue_cutoff = 8.41189e-06
# use the more precise one pvalues_ordered[20582] instead
# the difference is due to tiny computing precision issue
#  on the scale of e-20
# pvalues_ordered[20582] - pvalues_ordered[20583] = -1.863472e-20
pvalues_ordered = readRDS("../../data/intermediate_data/t8_outputs/pvalues_ordered.rds")
pvalue_cutoff = pvalues_ordered[20582]
asso_index_list = which(pvalues <= pvalue_cutoff)



# get the tcr name info

tcr_file = "t5_public_allele_level_tcr_filtered_wrt_vf_and_aa.csv"
public_allele_level_tcr_info = read.csv(paste0(data_dir, tcr_file), 
                                        header = TRUE)
n_tcr = dim(public_allele_level_tcr_info)[1]
n_tcr
#8739207

# get the hla name info

hla_file = "HLA_v2_features_reformat.csv"
HLA_v2_features = read.csv(paste0(data_dir, hla_file), 
                           header = TRUE)
n_hla = dim(HLA_v2_features)[1]
n_hla
# 215


# convert the asso index to the indexes of TCR and HLA

total_np = length(asso_index_list)
total_np
# 20582

tcr_list = rep("", total_np)
hla_allele_list = rep("", total_np)
association_pvalue_list  = rep(1.0, total_np)

for (k in 1:total_np){
  ind = asso_index_list[k]
  i = ceiling(ind/n_hla)
  j = ind - (n_hla * (i-1))
  cur_v_allele = public_allele_level_tcr_info$v_allele[i]
  cur_cdr3 = public_allele_level_tcr_info$amino_acids[i]
  cur_tcr = paste(cur_v_allele, ",", cur_cdr3, sep = "")
  cur_hla = HLA_v2_features$hla[j]
  if (pvalues[ind] != pvalues_big_list[[i]][[j]]){
    print(paste("indexes do not match at ind = ", ind, sep = ""))
    break
  } 
  tcr_list[k] = cur_tcr
  hla_allele_list[k] = cur_hla
  association_pvalue_list[k] = pvalues[ind]
}

# write out the asso pairs and p values
df_asso_v_allele = data.frame(tcr = tcr_list, 
                              hla_allele = hla_allele_list, 
                              association_pvalue = association_pvalue_list)

write.csv(df_asso_v_allele, 
          file = paste0(data_dir, "t12_HLA_associated_TCR_v_alleles.csv"), 
          row.names = FALSE)


sessionInfo()
q(save="no")

