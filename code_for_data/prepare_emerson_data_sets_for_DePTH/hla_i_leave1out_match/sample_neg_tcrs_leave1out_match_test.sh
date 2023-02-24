#!/bin/bash


for test_hla in HLA-B_08_01 HLA-B_07_02 HLA-C_07_01

do
    python sample_neg_tcrs_leave1out_match_test.py \
            --hla_class HLA_I \
            --test_hla $test_hla \
            --label leave1out_match \
            --n_fold 5
done
