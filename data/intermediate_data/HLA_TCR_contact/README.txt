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
