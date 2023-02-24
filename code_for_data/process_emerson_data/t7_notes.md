(0) V_allele level tcrs counting

file 't5_in_allele_level_tcr_dict.pickle' gives allele level tcr info for

69,122,667

unique tcrs.

File 't5_tcr_df_filter_out.py' filters out 878 tcrs by matching (v_family, amino_acid) tuple with those in file
"elife-38358-fig11-data1-v2.tds", leaving us

69,121,789

unique tcrs (saved as "t5_allele_level_tcr_filtered_wrt_vf_and_aa.csv"), including

8,739,207

unique public tcrs (saved as "t5_public_allele_level_tcr_filtered_wrt_vf_and_aa.csv").


(1) p-value computing


The new p value results contain 1,878,929,505 p values for 8,739,207 * 215 (v_allele_level_tcr, hla) pairs. The results are saved in rds format (list of lists, with the top layer of the list corresponding to tcrs, and the bottom layer of the list corresponding to hlas) (about 2GB) and saved in txt format (matrix format, with different rows corresponding to tcrs and different columns corresponding to hlas) (about 8GB). The tcr and hla names are not stored within the files due to the file sizes, but the tcr order follows that in "t5_public_allele_level_tcr_filtered_wrt_vf_and_aa.csv" and the hla order follows that in "HLA_v2_features_reformat.csv".
