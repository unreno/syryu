#!/usr/bin/env bash


script=`basename $0`
script_dir=`dirname $0`

#	Defaults:


function usage(){
#	echo
#	echo "Takes tab-separated file formatted like ..."
#	echo
#	echo "#SpecFile	SpecID	ScanNum	Title	FragMethod	Precursor	IsotopeError	PrecursorError(Da)	ChargePeptide	Protein	DeNovoScore	MSGFScore	SpecEValue	EValue	QValue	PepQValue"
#	echo "140521_EOC_MCis_T2_3.SCANS.mgf	index=5524	8449	140521_EOC_MCis_T2_3.8449.8449.2 File:\"140521_EOC_MCis_T2_3.raw\", NativeID:\"controllerType=0 controllerNumber=1 scan=8449\"	CID	1223.4911	0	0.0014648438	2	KVVDYSQFQES+79.966DDADEDYGR	sp|Q9H1E3|NUCKS_HUMAN(pre=R,post=D)	32	22	2.001993E-20	4.3754356E-13	0.0	0.0"
#	echo "140521_EOC_MCis_T2_3.SCANS.mgf	index=6007	9114	140521_EOC_MCis_T2_3.9114.9114.2 File:\"140521_EOC_MCis_T2_3.raw\", NativeID:\"controllerType=0 controllerNumber=1 scan=9114\"	CID	1159.4442	0	0.0020751953	2	VVDYS+79.966QFQESDDADEDYGR	sp|Q9H1E3|NUCKS_HUMAN(pre=K,post=D)	23	9	2.3498118E-16	5.1284723E-9	0.0	0.0"
#	echo "140521_EOC_MCis_T2_3.SCANS.mgf	index=6007	9114	140521_EOC_MCis_T2_3.9114.9114.2 File:\"140521_EOC_MCis_T2_3.raw\", NativeID:\"controllerType=0 controllerNumber=1 scan=9114\"	CID	1159.4442	0	0.0020751953	2	VVDY+79.966SQFQESDDADEDYGR	sp|Q9H1E3|NUCKS_HUMAN(pre=K,post=D)	23	9	2.3498118E-16	5.1284723E-9	0.0	0.0"
#	echo
#	echo "And adds 10 columns: 'mod1', 'mod2', ... , 'mod10'"
#	echo "These will contain the mass shift numbers for the modifications."
#	echo "Also adds Another column 'PlainPeptide', which contain the peptide sequence without modification."
#	echo "For example, 'SLDS+79.966DES+79.966EDEEDDYQQKR', mod1 has +79.966, mod2 has +79.966, and NA for the rest of mod columns."
#	echo "For 'PlainPeptide' column, it will contain 'SLDSDESEDEEDDYQQKR'."
	echo
	echo
	echo "1.    Read a peptide sequence and charge state from a formatted MODa’s output file (“out_blind_140521….txt.OUTPUT.tsv”)."
	echo "2.    For each sequence, calculate theoretical m/z values."
	echo "3.    For each entry in a formatted MODa’s output file, read the corresponding observed m/z values in the corresponding mgf file."
	echo "4.    For each sequence, compare theoretical m/z values and observed m/z values and calculate some score (I will give you score formula)."
	echo "5.    Attach the scores from #4 in a formatted MODa’s output file"
	echo ""
	echo "For #3-5, I will send you the detailed instruction after #1 & #2."
	echo "For now, I want the codes that do #1 and #2:"
	echo ""
	echo "For each line in a formatted MODa’s output file (“out_blind_140521….txt.OUTPUT.tsv”):"
	echo "A.   Read a peptide sequence and charge state."
	echo "(For example, the first line has peptide, “R.KTNDKDEKKEDGKQAENDSSNDDKTKK.S” with charge state 5 and the second line has “R.KTNDKDEKKEDGKQAENDSSNDDKT+128K.K” with a charge state 5.)"
	echo ""
	echo "B.   Using “theospec” codes (C codes), calculate theoretical m/z value for the sequence. “theospec” has a slightly different peptide sequence format than MODa output.  "
	echo "You need to either modify our peptide sequence in MODa or modify theospec codes."
	echo "What theospec can do:"
	echo "For the first sequence, I need a command “./theospec -z5 R.KTNDKDEKKEDGKQAENDSSNDDKTKK.S”"
	echo "For the second sequence, I need a command “./theospec -z5 -c*=128 R.KTNDKDEKKEDGKQAENDSSNDDKT*K.K”"
	echo "We may have more than one modification, so if we have “R.KTND-12KDEKKEDGKQAENDSSNDDKT+128K.K”,"
	echo "The command is like this: ./theospec -z1 -c*=128 -c#=-12 R.KTND#KDEKKEDGKQAENDSSNDDKT*K.K"
	echo ""
	echo "For now, printing output in the output file will do. But in the future, you will pass these output values to another function that you will write (#4 step)."
	echo "I attached out_blind_140521….txt.OUTPUT."
	echo "Theospec is in https://sourceforge.net/projects/protms/files/theospec/3.00/"
	echo
	echo
	echo "Usage:"
	echo
	echo "$script <OPTIONS> input_file(s)"
	echo
	echo "Will create each input file with the added extension .MSranker.tsv"
	echo
#	echo "Options:"
#	echo "	--verbose .................. NO VALUE, just boolean flag."
#	echo
#	echo "Default option values:"
#	echo "	--verbose ........... ${verbose}"
#	echo
#	echo "Notes:"
#	echo
#	echo "Examples:"
	echo
	echo
	exit 1
}



while [ $# -ne 0 ] ; do
	case $1 in
#		--min_shift)
#			shift; min_shift=$1; shift;;
#		--max_shift)
#			shift; max_shift=$1; shift;;
#		--min_shifted_mass)
#			shift; min_shifted_mass=$1; shift;;
#		--exclusion_min)
#			shift; exclusion_min=$1; shift;;
#		--exclusion_max)
#			shift; exclusion_max=$1; shift;;
#		--percent*)
#			shift; percent_of_peaks=$1; shift;;
#		-v|--v*)
#			verbose=true; shift;;
		-*)
			echo ; echo "Unexpected args from: ${*}"; usage ;;
		*)
			break;;
	esac
done

#	Basically, this is TRUE AND DO ...
[ $# -eq 0 ] && usage


while [ $# -ne 0 ] ; do

#	gawk -F"\t" '
#	MUST use gawk to use join
	perl -pe 's/\r//g' $1 | gawk -F"\t" '
	@include "join"
	BEGIN{ OFS="\t" }
	( NR == 1 ){ next }
	{
		split($11,parts,/\./);
		a=parts[1];
		b=join(parts,2,length(parts)-1,".");
		c=parts[length(parts)];

#		split($11,numbers,/[A-Z+]+/)
#		split($11,letters,/[0-9.+-]+/)
#		print $11, length(numbers), length(letters)
#		theospec_command="theospec -z"$4
#		print $11,$4,theospec_command
	}'
# > $1.peptide_and_charge.tsv
#' > $1.MSranker.tsv


#	These peptide sequences can start with a "-"??? Splitting on . and parse the middle.

#	I think that I could pipe the output from awk to bash to execute

#	For the final output file, appending MSranker.tsv will be good. For intermediate output files for step 1 and 2, filename_theoSpec.tsv will be good.






#' > $1.MSranker.tsv

	shift

done
