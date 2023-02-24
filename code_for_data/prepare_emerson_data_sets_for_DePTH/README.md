
# Prepare TCR-HLA pairs data sets from Emerson data for DePTH

There are three general settings, HLA-I, HLA-I leave1out, and HLA-II. The procedure of preparing Emerson (pos, neg) * (training, validation, test) files is similar for these settings, with more details for HLA-I leave1out in choosing TCR-HLA pairs.


[u1_cut_public_tcr_file.py] command line to cut public allele level TCR file

../../data/intermediate_data/t5_public_allele_level_tcr_filtered_wrt_vf_and_aa.csv to only keep the TCR names. The output file is used for all three settings.


[u2_v_gene_format_dict.ipynb] and manual procedure produce a file translating the v gene format from that in the extracted Emerson public TCRs to that in table

../../data/for_encoders/combo_xcr.tsv

The produced file is used by all HLA-I, HLA-I leave1out and HLA-II settings for translating the format of V gene.

Details for code files in each setting are put in the README.md of corresponding folder, with a more detailed version in HLA-I/HLA-II, and a shorter version in HLA-I leave1out.


<br />  
