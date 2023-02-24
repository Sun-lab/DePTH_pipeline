# this file is the third part for getting p-value cutoff for FDR 0.05

# for each of the 20 permutations, get the number of false positives at each true data 
# p value cutoff.
# first, sort the p values; O(nlogn)
# second, for each cutoff, do binary search. O(nlogn)

# this part is be done in the format of 20 job submissions, one for each permutation




args = commandArgs(trailingOnly = TRUE)
eval(parse(text=args[[1]]))


s = as.numeric(s)
s

# load ordered pvalues based on true data

sorted_file = "../../data/intermediate_data/t8_outputs/pvalues_ordered.rds"
pvalues_ordered = readRDS(sorted_file)

# load pvalues based on one permutation and convert into a vector
permutation_dir = "../../data/intermediate_data/t9_outputs/"
file_name = paste(permutation_dir, "pvalues_big_list_permuted_seed_", s, ".rds", sep = "")
per_pvalues_rds = readRDS(file_name)
per_pvalues = unlist(per_pvalues_rds)

rm(per_pvalues_rds)
gc()

per_pvalues_ordered = sort(per_pvalues)

rm(per_pvalues)
gc()


output_dir = "../../data/intermediate_data/t10_outputs/"
  
saveRDS(per_pvalues_ordered, 
        file = paste(output_dir, "permutation_pvalues_ordered_vec_seed_", s, ".rds", sep = ""))

length(pvalues_ordered) == length(per_pvalues_ordered)


# binary search to find the position of the last value <= the cutoff
cnt_leq <- function(per_pvalues_ordered, cutoff, i, j){
  if (per_pvalues_ordered[i] > cutoff){
    return(i-1)
  }
  if (per_pvalues_ordered[j] <= cutoff){
    return(j)
  }
  while (i< j){
    mid = floor((i+j)/2)
    if (per_pvalues_ordered[mid] <= cutoff){
      i = mid + 1
    }else if (per_pvalues_ordered[mid] > cutoff){
      j = mid
    }
  }
  return(i-1)
}


Sys.time()

cum_cnt_vec = rep(0, length(pvalues_ordered))
for (i in 1:length(pvalues_ordered)){
  cum_cnt_vec[i] = cnt_leq(per_pvalues_ordered, pvalues_ordered[i], 1, length(pvalues_ordered))
}

Sys.time()

saveRDS(cum_cnt_vec, 
        file = paste(output_dir, "cum_cnt_vec_seed_", s, ".rds", sep = ""))

sessionInfo()
q(save="no")

