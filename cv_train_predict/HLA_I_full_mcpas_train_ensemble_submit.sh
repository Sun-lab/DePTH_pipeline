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
  DePTH train --hla_class HLA_I \
              --data_dir data/HLA_I_full_mcpas/train_valid \
              --model_dir saved_models/HLA_I_full_mcpas_ensemble/HLA_I_full_mcpas_model_${rseed}_${np_seed}_${tf_seed} \
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
