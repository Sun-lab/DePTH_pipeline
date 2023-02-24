import os

cmd2 = "cut -f1,2 -d ',' ../../../data/intermediate_data/t5_public_allele_level_tcr_filtered_wrt_vf_and_aa.csv > ../../../data/intermediate_data/public_allele_level_tcr_name.csv"
os.system(cmd2)
