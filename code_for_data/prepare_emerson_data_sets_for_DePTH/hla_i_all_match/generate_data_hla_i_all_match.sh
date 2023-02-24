#!/bin/bash

python  generate_data.py \
        --hla_class HLA_I \
        --label all_match \
        --prop_train 0.6 \
        --prop_valid 0.2 \
        --n_fold 5 \
        --rseed 1627
