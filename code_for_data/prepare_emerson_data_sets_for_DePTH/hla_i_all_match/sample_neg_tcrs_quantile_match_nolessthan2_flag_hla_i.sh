#!/bin/bash



for i in {1..85}
do
  python sample_neg_tcrs_quantile_match_nolessthan2_flag.py \
          --hla_class HLA_I \
          --hla_i $i
done
