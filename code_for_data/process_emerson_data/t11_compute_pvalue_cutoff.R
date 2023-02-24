# this file is the fourth part for getting pvalue cutoff for FDR 0.05

# calculate 
#   \hat{q}(p_(n)) = E(# false positives with pvalue <= p_(n))/n; O(1)
# for i = n-1, n-2, ..., 1, calculate
#   \hat{q}(p_(i)) = 
#     min(E(# false positives with pvalue <= p_(i))/i, \hat{q}(p_(i+1)))); O(n)




# load the p values of the true data
sorted_file = "../../data/intermediate_data/t8_outputs/pvalues_ordered.rds"
pvalues_ordered = readRDS(sorted_file)
sum_cum_cnt_vec = rep(0, length(pvalues_ordered))



# load the number of false positives under all different pvalue cutoffs

rseed_vec = c(419, 581, 1242, 1590, 2175, 2282, 2834, 3539, 4567, 5270, 5532, 5560, 
              5561, 6022, 7035, 8249, 8303, 8439, 8581, 9510)

cnt_dir = "../../data/intermediate_data/t10_outputs/"

for (i in 1:20){
  s = rseed_vec[i]
  c4 = readRDS(paste(cnt_dir, "cum_cnt_vec_seed_", s, ".rds", sep = ""))
  print(s)
  print(length(c4))
  print(sum(c4 > 0))
  sum_cum_cnt_vec = sum_cum_cnt_vec + c4
}

average_cum_cnt_vec = sum_cum_cnt_vec/20

output_dir = "../../data/intermediate_data/t11_outputs/"
saveRDS(average_cum_cnt_vec, 
        file = paste0(output_dir, "average_cum_cnt_vec.rds"))


# compute estimated q values

n = length(pvalues_ordered)
q_hat_vec = rep(0, n)
q_hat_vec[n] = average_cum_cnt_vec[n]/n
for (j in 1:(n-1)){
  k = n - j
  q_hat_vec[k] = min(average_cum_cnt_vec[k]/k, q_hat_vec[k+1])
}

saveRDS(q_hat_vec, file = paste0(output_dir, "q_hat_vec.rds"))



# get the position of the last q value <= 0.05 cutoff

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

last_pos = cnt_leq(q_hat_vec, 0.05, 1, n)
# 20582
q_hat_vec[last_pos]
# 0.04978136
q_hat_vec[last_pos+1]
# 0.05012107


# get the corresponding p value cutoff
pvalue_cutoff = pvalues_ordered[last_pos]
# 8.41189e-06


