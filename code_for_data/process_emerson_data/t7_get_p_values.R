library(parallel)
library(doParallel)
library(doRNG)
library(foreach)


nCore = 4
#nCore = Sys.getenv("SLURM_CPUS_ON_NODE")

registerDoParallel(cores=nCore)
options(mc.cores=nCore)

RNGkind("L'Ecuyer-CMRG")

data_dir = "../../data/intermediate_data/"

tcr_file = "t5_public_allele_level_tcr_filtered_wrt_vf_and_aa.csv"
hla_file = "HLA_v2_features_reformat.csv"
  
public_allele_level_tcr_info = read.csv(paste0(data_dir, tcr_file), 
                                        header = TRUE)
dim(public_allele_level_tcr_info)[1]
head(public_allele_level_tcr_info)


HLA_v2_features = read.csv(paste0(data_dir, hla_file), 
                           header = TRUE)
dim(HLA_v2_features)


full_set = paste(c(0:665))

Sys.time()

pvalues_big_list = foreach(tcr_index = 1:dim(public_allele_level_tcr_info)[1]) %dorng%{
  
  tcr_p_set = strsplit(public_allele_level_tcr_info$individuals[tcr_index],",")[[1]]
  tcr_n_set = setdiff(full_set, tcr_p_set)
  
  pvalues_hla_slice = list()
  
  for (hla_index in 1:dim(HLA_v2_features)[1]){

    hla_name = HLA_v2_features$hla[hla_index]
    hla_p_set = strsplit(HLA_v2_features$ind_pos[hla_index], ",")[[1]]
    hla_n_set = strsplit(HLA_v2_features$ind_neg[hla_index], ",")[[1]]
    pp_count = length(intersect(tcr_p_set, hla_p_set))
    pn_count = length(intersect(tcr_p_set, hla_n_set))
    np_count = length(intersect(tcr_n_set, hla_p_set))
    nn_count = length(intersect(tcr_n_set, hla_n_set))
    
    testor = rbind(c(pp_count, pn_count), c(np_count, nn_count))
    pvalues_hla_slice[[hla_index]] = fisher.test(testor, alternative = 'greater')$p.value
  }
  pvalues_hla_slice
}

Sys.time()

saveRDS(pvalues_big_list, file = paste0(data_dir, "pvalues_big_list.rds"))

sessionInfo()
q(save="no")


