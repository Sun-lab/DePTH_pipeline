#!/bin/bash


for i in {1..85}
do
        R CMD BATCH "--args hla_i=${i}" tcr_freq_for_each_HLA_I.R tcr_freq_for_each_HLA_I_Rout/tcr_freq_for_each_HLA_I_${i}.Rout
done
