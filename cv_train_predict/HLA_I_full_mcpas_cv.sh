#!/bin/bash



declare -i i=0

for enc_method in one_hot blosum62 atchley pca; do
  for lr in 0.0001; do
    for n_dense in 1; do
      for n_units in [16] [32] [64]; do
        for dropout_flag in True; do
          for p_dropout in 0.2 0.5; do
            i+=1
            cur_name=${enc_method}_${n_dense}_${n_units}_${dropout_flag}_${p_dropout}
            DePTH cv --hla_class HLA_I \
                     --data_dir data/HLA_I_full_mcpas/train_valid \
                     --average_valid_dir results/cv_average_valid/HLA_I_full_mcpas/${cur_name} \
                     --enc_method ${enc_method} \
                     --lr ${lr} \
                     --n_dense ${n_dense} \
                     --n_units_str ${n_units} \
                     --dropout_flag ${dropout_flag} \
                     --p_dropout ${p_dropout} \
                     --rseed 1000 \
                     --np_seed 1216 \
                     --tf_seed 2207
          done
        done
      done
    done
  done
done


for enc_method in one_hot blosum62 atchley pca; do
  for lr in 0.0001; do
    for n_dense in 2; do
      for n_units in [32,16] [64,16]; do
        for dropout_flag in True; do
          for p_dropout in 0.2 0.5; do
            i+=1
            cur_name=${enc_method}_${n_dense}_${n_units}_${dropout_flag}_${p_dropout}
            DePTH cv --hla_class HLA_I \
                     --data_dir data/HLA_I_full_mcpas/train_valid \
                     --average_valid_dir results/cv_average_valid/HLA_I_full_mcpas/${cur_name} \
                     --enc_method ${enc_method} \
                     --lr ${lr} \
                     --n_dense ${n_dense} \
                     --n_units_str ${n_units} \
                     --dropout_flag ${dropout_flag} \
                     --p_dropout ${p_dropout} \
                     --rseed 1000 \
                     --np_seed 1216 \
                     --tf_seed 2207
          done
        done
      done
    done
  done
done


for enc_method in one_hot blosum62 atchley pca; do
  for lr in 0.0001; do
    for n_dense in 1; do
      for n_units in [16] [32] [64]; do
        for dropout_flag in False; do
          i+=1
          cur_name=${enc_method}_${n_dense}_${n_units}_${dropout_flag}
          DePTH cv --hla_class HLA_I \
                   --data_dir data/HLA_I_full_mcpas/train_valid \
                   --average_valid_dir results/cv_average_valid/HLA_I_full_mcpas/${cur_name} \
                   --enc_method ${enc_method} \
                   --lr ${lr} \
                   --n_dense ${n_dense} \
                   --n_units_str ${n_units} \
                   --dropout_flag ${dropout_flag} \
                   --rseed 1000 \
                   --np_seed 1216 \
                   --tf_seed 2207
        done
      done
    done
  done
done


for enc_method in one_hot blosum62 atchley pca; do
  for lr in 0.0001; do
    for n_dense in 2; do
      for n_units in [32,16] [64,16]; do
        for dropout_flag in False; do
          i+=1
          cur_name=${enc_method}_${n_dense}_${n_units}_${dropout_flag}
          DePTH cv --hla_class HLA_I \
                   --data_dir data/HLA_I_full_mcpas/train_valid \
                   --average_valid_dir results/cv_average_valid/HLA_I_full_mcpas/${cur_name} \
                   --enc_method ${enc_method} \
                   --lr ${lr} \
                   --n_dense ${n_dense} \
                   --n_units_str ${n_units} \
                   --dropout_flag ${dropout_flag} \
                   --rseed 1000 \
                   --np_seed 1216 \
                   --tf_seed 2207
        done
      done
    done
  done
done
