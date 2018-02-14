#!/usr/bin/env bash



#	I need a software (Modification Locator) that locates peptide modification in protein.
#
#	For evidence file, please use “Modified sequence” column and “Leading razor protein” column.


#==> evidence.txt <==
#Sequence	Length	Modifications	Modified sequence	Oxidation (M) Probabilities	Phospho (STY) Probabilities	Oxidation (M) Score Diffs	Phospho (STY) Score Diffs	Acetyl (Protein N-term)	Oxidation (M)	Phospho (STY)	Missed cleavages	Proteins	Leading proteins	Leading razor protein	Type	Raw file	Fraction	Experiment	MS/MS m/z	Charge	m/z	Mass	Resolution	Uncalibrated - Calibrated m/z [ppm]	Uncalibrated - Calibrated m/z [Da]	Mass Error [ppm]	Mass Error [Da]	Uncalibrated Mass Error [ppm]	Uncalibrated Mass Error [Da]	Max intensity m/z 0	Retention time	Retention length	Calibrated retention time	Calibrated retention time start	Calibrated retention time finish	Retention time calibration	Match time difference	Match m/z difference	Match q-value	Match score	Number of data points	Number of scans	Number of isotopic peaks	PIF	Fraction of total spectrum	Base peak fraction	PEP	MS/MS Count	MS/MS Scan Number	Score	Delta score	Combinatorics	Intensity	Reverse	Potential contaminant	id	Protein group IDs	Peptide ID	Mod. peptide ID	MS/MS IDs	Best MS/MS	AIF MS/MS IDs	Oxidation (M) site IDs	Phospho (STY) site IDs

#	head evidence.txt | awk -F"\t" '{print $4,$15}'
#Modified sequence Leading razor protein
#_AAAAAAAAAAAAAAAGAGAGAK_ tr|Q53ZR1|Q53ZR1_HUMAN
#_AAAAAAAAAAAAAAAGAGAGAK_ tr|Q53ZR1|Q53ZR1_HUMAN
#_AAAAAAAAAAAAAAAGAGAGAK_ tr|Q53ZR1|Q53ZR1_HUMAN
#_AAAAAAAAAAAAAAAGAGAGAK_ tr|Q53ZR1|Q53ZR1_HUMAN
#_AAAAAAAAAAAATGTEAGPGTAGGSENGSEVAAQPAGLSGPAEVGPGAVGER_ sp|O60341|KDM1A_HUMAN
#_AAAAAAAAAAAATGTEAGPGTAGGSENGSEVAAQPAGLSGPAEVGPGAVGER_ sp|O60341|KDM1A_HUMAN
#_AAAAAAAAAAAATGTEAGPGTAGGSENGSEVAAQPAGLSGPAEVGPGAVGER_ sp|O60341|KDM1A_HUMAN
#_AAAAAAAAAAAATGTEAGPGTAGGSENGSEVAAQPAGLSGPAEVGPGAVGER_ sp|O60341|KDM1A_HUMAN
#_AAAAAAAAAAAATGTEAGPGTAGGSENGSEVAAQPAGLSGPAEVGPGAVGER_ sp|O60341|KDM1A_HUMAN

#
#	I uploaded “evidence” file and “PeptideMatach” software in Azure Modification Locator folder. I uploaded fasta file which contains protein sequences in “Blob container”. (I failed to upload it to the same folder.)
#
#	 
#
#	I want two output files contains:
#
#	The first output file (name it as MatchedModification.txt) contains 1) Modified sequence, 2) Leading razor protein, 3) Start position, 4) End position, 5) STY position, 6) STY modification position, 7) M position, and 8) M modification position. Please give users option that they can choose their amino acids of interest. In this example, users chose “STY”and “M”.
#	For the given “Modified sequence”, please fine the corresponding “Leading razor protein” in fasta file.
#
#	Then, find their start/end position of sequence in protein. The first position of protein starts with 1 instead of 0.
#
#	(I uploaded that codes (PeptideMatch java codes) that does this. It was previously published by another group. If you like, you can modify these codes or re-write it in your language you prefer.)
#
#	                If “Modified sequence” contains S, T, or Y, then please look for the location of either S, T, or Y in protein and record it. If there are more than one of these, then record all separated by comma.
#
#	                If S, T, or Y are modified in “Modified sequence”, then please record its location in STY modification position column.
#
#	                If “Modified sequence” contains M, then please look for the location of M in protein and record it. If there are more than one site, then record them separated by comma.
#
#	                If M modified in “Modified sequence”, then please record its location in M modification position column.
#
#	 
#
#	The second output file (name it as ProteinModification.txt) contains 1) Leading razor protein, 2) AA, 3) AA location, 4) Sequences with Modified AA, and 5) Sequences with Unmodified AA
#	For each leading razor protein that were observed in at least once in evidence file, find all S, T, Y, or M. AA will contain one of these amino acids located in a specific location.
#
#	Each location of S, T, Y, or M will be one line. “AA location” will contain the location of the corresponding amino acid (one location per line).
#
#	For “Sequences with Modified AA”, records all unique sequences with the corresponding modified “AA”.
#
#	For “Sequences with Unmodified AA”, records all unique sequences with the corresponding unmodified “AA”.
#
#	 
#
#	 
#
#	Cartoon Sample
#
#	Fasta file:
#
#	>PROTEIN1
#
#	AAASSSTTYAAAM
#
#	 
#
#	Evidence file:
#
#	Modified Sequence Leading razor protein
#
#	_ASSST[80]TY_
#
#	_SSS[80]T[80]_
#
#	 
#
#	Then,
#
#	Output file 1:
#
#	Modified sequence\t Leading razor protein\t Start position\t End position\t STY position\t STY modification position\t M position\t M modification position \n
#
#	_ASSST[80]TY_\t       PROTEIN1\t     3\t      9\t    4;5;6;7;8;9\t    7\t        \t      \t
#
#	_SSS[80]T[80]_\t      PROTEIN1\t     4\t      7\t    4;5;6;7\t        4;7\t      \t      \t
#
#	 
#
#	Output file 2:
#
#	Leading razor protein\t AA\t AA location\t Sequences with Modified AA\t Sequences with Unmodified AA\t
#
#	PROTEIN1\t       S\t    4\t     \t      _ASSST[80]TY_;_SSS[80]T[80]_
#
#	PROTEIN1\t       S\t    5\t     \t      _ASSST[80]TY_;_SSS[80]T[80]_
#
#	PROTEIN1\t       S\t    6\t     \t      _SSS[80]T[80]_\t                 _ASSST[80]TY_\t
#
#	PROTEIN1\t       T\t    7\t\    \t      _ASSST[80]TY_; _SSS[80]T[80]_\t                   \t                                                            
#
#	…
#
#	 
#
#	Thank you!!
#
#
#	I'll see if I can just uniq the modified sequence / leading razor protein from within the script itself.
#
#	All positions are relative to the start of the full protein found in the fasta file.
#
#	Ignore the underscores.
#
#	Modifications come after the column.
#
#	The modification value itself is irrelevant in this script.




#	Select only those matching the requested amino acids.

#	[[ -f evidence.STYM.txt ]] || awk -F"\t" '( $4 ~ /[STYM]/ ){print $4,$15}' evidence.txt | sort -u > evidence.STYM.txt


#	This will likely be better done in something like ruby.





