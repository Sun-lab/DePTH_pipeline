#!/bin/bash



test_hla_i=1

for hla_i in {2..85}
do
  python sample_neg_tcrs_leave1out_match_train_valid.py \
            --hla_class HLA_I \
            --test_hla_i ${test_hla_i} \
            --hla_i ${hla_i} \
            --label leave1out_match \
            --n_fold 5
done


test_hla_i=11

for hla_i in {1..10} {12..85}
do
  python sample_neg_tcrs_leave1out_match_train_valid.py \
            --hla_class HLA_I \
            --test_hla_i ${test_hla_i} \
            --hla_i ${hla_i} \
            --label leave1out_match \
            --n_fold 5
done


test_hla_i=23

for hla_i in {1..22} {24..85}
do
  python sample_neg_tcrs_leave1out_match_train_valid.py \
            --hla_class HLA_I \
            --test_hla_i ${test_hla_i} \
            --hla_i ${hla_i} \
            --label leave1out_match \
            --n_fold 5
done
