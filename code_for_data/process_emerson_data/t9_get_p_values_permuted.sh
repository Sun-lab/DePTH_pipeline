#!/bin/bash

for seed in $( cat seeds_20.txt)
do
    R CMD BATCH '--args s='$seed'' t9_get_p_values_permuted.R t9_get_p_values_permuted_$seed.Rout
done
