#!/bin/bash

for seed in $( cat seeds_20.txt)
do
    R CMD BATCH '--args s='$seed'' t10_get_counts.R t10_get_counts_$seed.Rout
done
