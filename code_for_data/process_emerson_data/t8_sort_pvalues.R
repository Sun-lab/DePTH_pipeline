# this file is the first part for getting pvalue cutoff for FDR 0.05

# sort the p values based on true data; O(nlogn) 

# load pvalues based on true data

data_dir = "../../data/intermediate_data/"


pvalues_rds = readRDS(paste0(data_dir, "pvalues_big_list.rds"))
pvalues = unlist(pvalues_rds)
pvalues_ordered = sort(pvalues)

output_dir = "../../data/intermediate_data/t8_outputs/"
saveRDS(pvalues_ordered, file = paste0(output_dir, "pvalues_ordered.rds"))













