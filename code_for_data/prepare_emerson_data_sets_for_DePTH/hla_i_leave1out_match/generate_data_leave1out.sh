#!/bin/bash


for test_hla_i in 1 11 23
do
  python generate_data_leave1out.py \
            --hla_class HLA_I \
            --test_hla_i ${test_hla_i} \
            --limit_hla_i 85 \
            --label leave1out_match \
            --prop_train 0.75 \
            --n_fold 5 \
            --rseed 1627
done
