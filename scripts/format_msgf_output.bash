#!/usr/bin/env bash


script=`basename $0`
script_dir=`dirname $0`

#	Defaults:




function usage(){
	echo
	echo "Takes tab-separated file formatted like ..."
	echo
	echo "#SpecFile	SpecID	ScanNum	Title	FragMethod	Precursor	IsotopeError	PrecursorError(Da)	ChargePeptide	Protein	DeNovoScore	MSGFScore	SpecEValue	EValue	QValue	PepQValue"
	echo "140521_EOC_MCis_T2_3.SCANS.mgf	index=5524	8449	140521_EOC_MCis_T2_3.8449.8449.2 File:\"140521_EOC_MCis_T2_3.raw\", NativeID:\"controllerType=0 controllerNumber=1 scan=8449\"	CID	1223.4911	0	0.0014648438	2	KVVDYSQFQES+79.966DDADEDYGR	sp|Q9H1E3|NUCKS_HUMAN(pre=R,post=D)	32	22	2.001993E-20	4.3754356E-13	0.0	0.0"
	echo "140521_EOC_MCis_T2_3.SCANS.mgf	index=6007	9114	140521_EOC_MCis_T2_3.9114.9114.2 File:\"140521_EOC_MCis_T2_3.raw\", NativeID:\"controllerType=0 controllerNumber=1 scan=9114\"	CID	1159.4442	0	0.0020751953	2	VVDYS+79.966QFQESDDADEDYGR	sp|Q9H1E3|NUCKS_HUMAN(pre=K,post=D)	23	9	2.3498118E-16	5.1284723E-9	0.0	0.0"
	echo "140521_EOC_MCis_T2_3.SCANS.mgf	index=6007	9114	140521_EOC_MCis_T2_3.9114.9114.2 File:\"140521_EOC_MCis_T2_3.raw\", NativeID:\"controllerType=0 controllerNumber=1 scan=9114\"	CID	1159.4442	0	0.0020751953	2	VVDY+79.966SQFQESDDADEDYGR	sp|Q9H1E3|NUCKS_HUMAN(pre=K,post=D)	23	9	2.3498118E-16	5.1284723E-9	0.0	0.0"
	echo
	echo "And adds 10 columns: 'mod1', 'mod2', ... , 'mod10'"
	echo "These will contain the mass shift numbers for the modifications."
	echo "Also adds Another column 'PlainPeptide', which contain the peptide sequence without modification."
	echo "For example, 'SLDS+79.966DES+79.966EDEEDDYQQKR', mod1 has +79.966, mod2 has +79.966, and NA for the rest of mod columns."
	echo "For 'PlainPeptide' column, it will contain 'SLDSDESEDEEDDYQQKR'."
	echo
	echo
	echo "Usage:"
	echo
	echo "$script <OPTIONS> input_file(s)"
	echo
	echo "Will create each input file with the added extension .OUTPUT.tsv"
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



# Modify your “format_moda_output.bash” to add one more column (i.e. “PlainPeptide”)  to the output file. For example, if the peptide sequence is “K.T+128NDKDEKKE.S”, then “PlainPeptide” column will have any sequence between “.” w/o modification information. Thus, the it will contain “TNDKDEKKE”. The example input file is attached (out_param_140521_EOC_MCis_T2_3.SCANS.mgf.txt).



while [ $# -ne 0 ] ; do

#	Remove possible Windows ^Ms
#	doesn't work on older versions of sed.
#sed 's/\r/\n/g' $1 | awk '

#perl -pe 's/\r/\n/g' $1 | awk '


#	NOTE THAT THIS TIME THE mods can INCLUDE DECIMAL POINTS

#	why not have -F"\t"?

perl -pe 's/\r//g' $1 | awk '
BEGIN {
	FS="\t";
	OFS="\t";
	mods=10;
}
( NR == 1 ){
	out=$0
	for(i=1;i<=mods;i++){ out=out OFS "modAA"i OFS "mod"i }
	print out,"PlainPeptide"
	next
}
{
	out=$0

	split($10,a,/[A-Z+]*/)
	split($10,b,/[0-9+-.]*/)
	for(i=1;i<=mods;i++){
		modAA=(length(b) > i) ? substr(b[i],length(b[i]),1) : "NA"
		mod=(length(a) > i+1) ? a[i+1] : "NA"
		out=out OFS modAA OFS mod
	}

	PlainPeptide=$10
	gsub(/[0-9+-.]*/,"",PlainPeptide)
	out=out OFS PlainPeptide
	
	print out
}
' > $1.WITH_MODS.tsv

	shift

done
