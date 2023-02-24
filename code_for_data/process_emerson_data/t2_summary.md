This file summarizes what jupyter notebook t2_check_additional_pos_contacts.log.ipynb does and finds.

"t2_check_additional_pos_contacts.log.ipynb" explores the contact position information  from the file "../../data/intermediate_data/HLA_TCR_contact/contacts.log", in order to find possible additional positions to add to the 34 peptide contact positions used for HLA-I alleles, the 15 positions used for HLA-II alpha chain, and the 19 positions for HLA-II beta chains.

The 34 positions (1-indexed) for HLA-I alleles from NetMHCPan are:

[7, 9, 24, 45, 59, 62, 63, 66, 67, 69, 70, 73, 74, 76, 77, 80, 81, 84, 95, 97, 99, 114, 116, 118, 143, 147, 150, 152, 156, 158, 159, 163, 167, 171]

The 15 positions (1-indexed, but needs other additional adjustment) for HLA-II alleles from NetMHCIIpan-3.0 are:

[9, 11, 22, 24, 31, 52, 53, 58, 59, 61, 65, 66, 68, 72, 73]

The 19 positions (1-indexed, but needs other additional adjustment) for HLA-II alleles from NetMHCIIpan-3.0 are:

[9, 11, 13, 26, 28, 30, 47, 57, 67, 70, 71, 74, 77, 78, 81, 85, 86, 89, 90]


The contacts table is separated into six data frames:

tcr_HLA_I_contacts

tcr_HLA_II_contacts_alpha

tcr_HLA_II_contacts_beta

pep_HLA_I_contacts

pep_HLA_II_contacts_alpha

pep_HLA_II_contacts_beta

A table of contact positions and corresponding counts is written out for each data frame. The six output files are:

"../../data/intermediate_data/step74_tcr_HLA_I_contacts.csv"

"../../data/intermediate_data/step74_pep_HLA_I_contacts.csv"

"../../data/intermediate_data/step74_tcr_HLA_II_contacts_alpha.csv"

"../../data/intermediate_data/step74_pep_HLA_II_contacts_alpha.csv"

".../../data/intermediate_data/step74_tcr_HLA_II_contacts_beta.csv"

"../../data/intermediate_data/step74_pep_HLA_II_contacts_beta.csv"

(these files are not used in later steps for adding positions, as later steps directly copy and paste the selected positions from the notebook output.)

For tcr_HLA_I_contacts, we take a subset of all contact positions with count>=10, and do the same for pep_HLA_I_contacts. Then we union these two subsets together, and get the positions that fall into this union but not in the original 34 positions (after changed to 0-indexed). This give us 17 additional positions: [4, 57, 64, 67, 71, 74, 122, 145, 148, 150, 153, 154, 156, 160, 161, 165, 169].

For tcr_HLA_II_contacts_alpha, we take a subset of all contact positions with count>=3, and do the same for pep_HLA_II_contacts_alpha (the cutoff 3 was chosen based on the min count of the 15 positions from NetMHCIIpan-3.0--after subtracting the positions by 4 all together--in the counter dictionary corresponding to pep_HLA_II_contacts_alpha). Then we union these two subsets together, and get the positions that fall into this union but not in the original 15 positions (after subtracting all 15 positions by 4 together). This gives us 18 additional positions: [4, 28, 35, 39, 45, 46, 47, 50, 51, 53, 56, 58, 60, 63, 65, 67, 71, 72].

For tcr_HLA_II_contacts_beta, we take a subset of all contact positions with count>=10, and do the same for pep_HLA_II_contacts_beta (the cutoff 10 was chosen based on the min count of the 19 positions from NetMHCIIpan-3.0--after subtracting the positions by 7 all together--in the counter dictionary corresponding to pep_HLA_II_contacts_alpha). Then we union these two subsets together, and get the positions that fall into this union but not in the original 19 positions (after subtracting all 19 positions by 7 together). This gives us 9 additional positions: [49, 53, 54, 57, 59, 62, 66, 69, 75].
