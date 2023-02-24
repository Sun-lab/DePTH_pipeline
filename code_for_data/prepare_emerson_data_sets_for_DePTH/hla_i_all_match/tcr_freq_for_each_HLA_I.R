
# For each HLA, 
# this code outputs the TCRs that appear at least twice among the subjects 
# with known HLA existing or not status
# and their corresponding frequency among the subjects under consideration

# the public TCR and HLA allele information used here are based on all 666 
# Emerson HIP subjects

args = commandArgs(trailingOnly=TRUE)
args

if (length(args) != 1) {
  message("one argument is expected, use 1 as default.\n")
  hla_i = 1
}else{
  eval(parse(text=args[[1]]))
}

hla_i = as.integer(hla_i)


data_dir = "../../../data/intermediate_data/"

public_allele_level_tcr_info = 
  read.csv(paste0(data_dir, 
                  "t5_public_allele_level_tcr_filtered_wrt_vf_and_aa.csv"), 
           header = TRUE)
dim(public_allele_level_tcr_info)[1]
head(public_allele_level_tcr_info)



HLA_v2_features = 
  read.csv(paste0(data_dir, "HLA_v2_features_reformat.csv"), header = TRUE)
dim(HLA_v2_features)


# pick out the 85 HLA-I alleles

hla_short = rep(NA, dim(HLA_v2_features)[1])

for (i in 1:dim(HLA_v2_features)[1]){
  hla_short[i] = substr(HLA_v2_features$hla[i], 1, 5)
}

HLA_v2_features$hla_short = hla_short

HLA_I_v2_features = 
  HLA_v2_features[HLA_v2_features$hla_short %in% c("HLA-A", "HLA-B", "HLA-C"),]

dim(HLA_I_v2_features)


full_set = paste(c(0:665))


  
hla_p_set = strsplit(HLA_I_v2_features$ind_pos[hla_i], ",")[[1]]
hla_n_set = strsplit(HLA_I_v2_features$ind_neg[hla_i], ",")[[1]]

pp_cnt = rep(NA, dim(public_allele_level_tcr_info)[1])
pn_cnt = rep(NA, dim(public_allele_level_tcr_info)[1])
np_cnt = rep(NA, dim(public_allele_level_tcr_info)[1])
nn_cnt = rep(NA, dim(public_allele_level_tcr_info)[1])

for (tcr_i in 1:dim(public_allele_level_tcr_info)[1]){
  
  tcr_p_set = strsplit(public_allele_level_tcr_info$individuals[tcr_i],",")[[1]]
  tcr_n_set = setdiff(full_set, tcr_p_set)

  pp_cnt[tcr_i] = length(intersect(tcr_p_set, hla_p_set))
  pn_cnt[tcr_i] = length(intersect(tcr_p_set, hla_n_set))
  np_cnt[tcr_i] = length(intersect(tcr_n_set, hla_p_set))
  nn_cnt[tcr_i] = length(intersect(tcr_n_set, hla_n_set))
  
}

p_cnt = pp_cnt + pn_cnt

chosen_tcr_i = which(p_cnt>=2)

p_cnt_chosen = p_cnt[chosen_tcr_i]
v_allele_chosen = public_allele_level_tcr_info$v_allele[chosen_tcr_i]
aas_chosen = public_allele_level_tcr_info$amino_acids[chosen_tcr_i]  

df_cur = data.frame(v_allele = v_allele_chosen, 
                    amino_acids = aas_chosen, 
                    freq = p_cnt_chosen)

fname = paste0("HLA_I_all_match_prepare/tcr_freq_HLA_I_",
               as.character(hla_i), ".csv")

write.csv(df_cur, 
	  file = paste0(data_dir, fname), 
	  row.names = FALSE, 
	  quote=FALSE)






gc()

sessionInfo()
q(save="no")
