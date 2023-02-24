#!/bin/bash

for i in {1..125}
do
  python sample_neg_tcrs_quantile_match.py \
          --hla_class HLA_II \
          --hla_i $i
done
