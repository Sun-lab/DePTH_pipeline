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

  DePTH predict --test_file data/HLA_II_all_match/test/test_pos.csv \
                --hla_class HLA_II \
                --output_dir results/predicted_scores/HLA_II_all_match_ensemble/HLA_II_all_match_pos_test_${rseed}_${np_seed}_${tf_seed} \
                --default_model False \
                --model_dir saved_models/HLA_II_all_match_ensemble/HLA_II_all_match_model_${rseed}_${np_seed}_${tf_seed} \
                --enc_method one_hot

  DePTH predict --test_file data/HLA_II_all_match/test/test_neg.csv \
                --hla_class HLA_II \
                --output_dir results/predicted_scores/HLA_II_all_match_ensemble/HLA_II_all_match_neg_test_${rseed}_${np_seed}_${tf_seed} \
                --default_model False \
                --model_dir saved_models/HLA_II_all_match_ensemble/HLA_II_all_match_model_${rseed}_${np_seed}_${tf_seed} \
                --enc_method one_hot


done < "$seed_file"
