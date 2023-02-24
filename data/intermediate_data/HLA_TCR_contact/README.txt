Files in this folder are the pMHC-TCR contacts information from Phil Bradley. 

These results were derived from parsing TCR:pMHC structures and assigning contacts to the HLA positions. The tricky part is aligning the sequence in the PDB structure files to the HLA sequences in a consistent way. 

The files in this folder have the alignments for the different alleles and a log file with all the contacts. The log file contains lots of lines like the following two, with the columns:

1: tag should start either tcr or pep
2. pdb chain for mhc
3. pdb chain for other
4. allele
5. HLA alignment column making the contact, 0-indexed, numbered with respect to the alignments in the zip folder.
6. HLA sequence at this position
7. all HLA sequences in this column (you could use this to sanity-check the alignment)
8. all HLA sequences in this column seen in the Emerson dataset
9. PDB filename 

tcr_contact_pos: A D A*02 169 R IR R /home/pbradley/tcr_scripts/pdb_files/1ao7.pdb.human.MH1.A-02.A.C.DE.pdb
pep_contact_pos: A C A*02 4 M LMV M /home/pbradley/tcr_scripts/pdb_files/1ao7.pdb.human.MH1.A-02.A.C.DE.pdb

tcr_contact_pos is where HLA has contact with TCR and pep_contact_pos is where HLA has contact with peptide


It should also be possible to report on the TCR/peptide position numbers that are making these contacts to HLA. There again the tricky part is aligning from the structure files (where the numbering is pretty arbitrary) to a standard reference frame.

Q/A:

1. What is pdb chain for other, which are labeled by D and C in your two examples?

it's the chain assigned to the peptide (often C) or tcr (often D or E) in the pdb file. But it's not standardized at all: in some pdb files it might be A/B/H/K... 


2. The column “7. all HLA sequences in this column” is all the possible amino acids at that position for that HLA allele, e.g., A*02 in the two examples?

Yes, that's right. So column 7 includes all the crazy alleles like A*02:01:123 and A*02:171 whereas column 8 is just the A*02 ones that were seen in the Emerson CMV cohort of 666 individuals, so probably mainly A*02:01... 
 