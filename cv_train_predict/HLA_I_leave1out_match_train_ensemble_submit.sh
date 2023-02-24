#!/bin/bash

declare -i i=0


seed_file="ensemble_seeds_20.txt"
while IFS= read -r line
do
  stringarray=($line)
  rseed=${stringarray[0]}
  np_seed=${stringarray[1]}
  tf_seed=${stringarray[2]}

  i+=1

  leftout=B0702
  DePTH train --hla_class HLA_I \
              --data_dir data/HLA_I_leave1out_match/HLA_I_leave1out_match_${leftout}/train_valid \
              --model_dir saved_models/HLA_I_leave1out_match_${leftout}_ensemble/HLA_I_leave1out_match_${leftout}_model_${rseed}_${np_seed}_${tf_seed} \
              --enc_method pca \
              --lr 0.0001 \
              --n_dense 1 \
              --n_units_str [64] \
              --dropout_flag True \
              --p_dropout 0.5 \
              --rseed ${rseed} \
              --np_seed ${np_seed} \
              --tf_seed ${tf_seed}

    leftout=B0801
    DePTH train --hla_class HLA_I \
                --data_dir data/HLA_I_leave1out_match/HLA_I_leave1out_match_${leftout}/train_valid \
                --model_dir saved_models/HLA_I_leave1out_match_${leftout}_ensemble/HLA_I_leave1out_match_${leftout}_model_${rseed}_${np_seed}_${tf_seed} \
                --enc_method pca \
                --lr 0.0001 \
                --n_dense 1 \
                --n_units_str [64] \
                --dropout_flag True \
                --p_dropout 0.5 \
                --rseed ${rseed} \
                --np_seed ${np_seed} \
                --tf_seed ${tf_seed}

    leftout=C0701
    DePTH train --hla_class HLA_I \
                --data_dir data/HLA_I_leave1out_match/HLA_I_leave1out_match_${leftout}/train_valid \
                --model_dir saved_models/HLA_I_leave1out_match_${leftout}_ensemble/HLA_I_leave1out_match_${leftout}_model_${rseed}_${np_seed}_${tf_seed} \
                --enc_method pca \
                --lr 0.0001 \
                --n_dense 2 \
                --n_units_str [64,16] \
                --dropout_flag True \
                --p_dropout 0.2 \
                --rseed ${rseed} \
                --np_seed ${np_seed} \
                --tf_seed ${tf_seed}


done < "$seed_file"
