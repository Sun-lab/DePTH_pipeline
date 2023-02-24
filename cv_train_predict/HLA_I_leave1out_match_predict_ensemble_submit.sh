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

  for leftout in B0702 B0801 C0701; do

  DePTH predict --test_file data/HLA_I_leave1out_match/HLA_I_leave1out_match_${leftout}/test/test_pos.csv \
                --hla_class HLA_I \
                --output_dir results/predicted_scores/HLA_I_leave1out_match_${leftout}_ensemble/HLA_I_leave1out_match_${leftout}_pos_test_${rseed}_${np_seed}_${tf_seed} \
                --default_model False \
                --model_dir saved_models/HLA_I_leave1out_match_${leftout}_ensemble/HLA_I_leave1out_match_${leftout}_model_${rseed}_${np_seed}_${tf_seed} \
                --enc_method pca

  DePTH predict --test_file data/HLA_I_leave1out_match/HLA_I_leave1out_match_${leftout}/test/test_neg.csv \
                --hla_class HLA_I \
                --output_dir results/predicted_scores/HLA_I_leave1out_match_${leftout}_ensemble/HLA_I_leave1out_match_${leftout}_neg_test_${rseed}_${np_seed}_${tf_seed} \
                --default_model False \
                --model_dir saved_models/HLA_I_leave1out_match_${leftout}_ensemble/HLA_I_leave1out_match_${leftout}_model_${rseed}_${np_seed}_${tf_seed} \
                --enc_method pca

  done
done < "$seed_file"
