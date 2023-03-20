#!/bin/bash


seed_file="ensemble_seeds.txt"

while IFS= read -r line
do
  stringarray=($line)
  rseed=${stringarray[0]}
  np_seed=${stringarray[1]}
  tf_seed=${stringarray[2]}


  DePTH predict --test_file data/HLA_I_all_match/train_valid/valid_pos.csv \
                --hla_class HLA_I \
                --output_dir results/predicted_scores/HLA_I_all_match_valid_ensemble/HLA_I_all_match_pos_valid_${rseed}_${np_seed}_${tf_seed} \
                --default_model False \
                --model_dir saved_models/HLA_I_all_match_ensemble/HLA_I_all_match_model_${rseed}_${np_seed}_${tf_seed} \
                --enc_method one_hot

  DePTH predict --test_file data/HLA_I_all_match/train_valid/valid_neg.csv \
                --hla_class HLA_I \
                --output_dir results/predicted_scores/HLA_I_all_match_valid_ensemble/HLA_I_all_match_neg_valid_${rseed}_${np_seed}_${tf_seed} \
                --default_model False \
                --model_dir saved_models/HLA_I_all_match_ensemble/HLA_I_all_match_model_${rseed}_${np_seed}_${tf_seed} \
                --enc_method one_hot


done < "$seed_file"
