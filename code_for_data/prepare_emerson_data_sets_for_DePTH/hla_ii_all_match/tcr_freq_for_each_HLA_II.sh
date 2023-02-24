#!/bin/bash


for i in {1..125}
do
  R CMD BATCH "--args hla_i=${i}" tcr_freq_for_each_HLA_II.R tcr_freq_for_each_HLA_II_Rout/tcr_freq_for_each_HLA_II_${i}.Rout
done
