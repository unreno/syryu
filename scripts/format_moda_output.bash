#!/usr/bin/env bash


script=`basename $0`
script_dir=`dirname $0`

#	Defaults:





function usage(){
	echo
	echo "Takes tab-separated file formatted like ..."
	echo
	echo ">>K:/Phospho_Goodlett/raw/140521_EOC_MCis_T2_3.mgf	14	2869.4179	6	0		"
	echo "2870.6817	-1.2638	22	0.0004	R.A+5RVILGSPRPRVIVSSPWPAVVVASPR.P	sp|P0DI83	72~98"
	echo "2870.4205	-1.0026	18	0	R.NQCPAKPDGGGAPNGVRNGLAAELGPASPR.R	sp|Q9BX95	74~103"
	echo "2868.3922	1.0257	15	0	R.QEPALRGSPGPLTPHPCNELGP+94PASPR.T	sp|Q8TAY7	39~65"
	echo "2869.4381	-0.0202	14	0	K.KKGPPPIR-5SCDDFSHMGTLPHSKSPR.Q	sp|O75815	55~80"
	echo "2870.4587	-1.0408	13	0	R.PKAPASP-37RRPQTPTPSEQDADPGPASPR.D	sp|Q8NBB4	5~32"
	echo
	echo ">>K:/Phospho_Goodlett/raw/140521_EOC_MCis_T2_3.mgf	19	3094.4649	5	0		"
	echo "3094.465	-0.0001	82	0.9997	R.KTNDKDEKKEDGKQAENDSSNDDKTKK.S	sp|Q9BXP5	301~327"
	echo "3094.3701	0.0949	75	0.9971	R.KTNDKDEKKEDGKQAENDSSNDDKT+128K.K	sp|Q9BXP5	301~326"
	echo "3094.3701	0.0949	71	0.9929	K.T+128NDKDEKKEDGKQAENDSSNDDKTKK.S	sp|Q9BXP5	302~327"
	echo "3094.5495	-0.0845	46	0.0083	K.VLL-33AWSGGPSSS+55SMVWQVLEGLSQDSAKR-2.L	sp|Q2VPK5	75~103"
	echo "3093.4716	0.9933	45	0.0137	K.KTL-41SQMSLSSDNSHATQNISP+78KKDD+50FK.N	sp|Q6ZP01	731~757"
	echo
	echo "Each group starts with a line containing a filename in column 1 (ignored)"
	echo "	column 2 : spec_index"
	echo "	column 3 : observed_MW"
	echo "	column 4 : charge_state"
	echo "	column 5 : scan_number"
	echo "Each line after, but before a blank line, is in order of rank and contains ..."
	echo "	column 1 : calculated_MW"
	echo "	column 2 : delta_mass"
	echo "	column 3 : score"
	echo "	column 4 : probability"
	echo "	column 5 : peptide"
	echo "	column 6 : protein"
	echo "	column 7 : pept_position"
	echo
	echo "The peptide string can contain multiple (up to 3) + or - INTEGERs to be selected as mod1,mod2 and mod3."
	echo "These should be NA when not set."
	echo
	echo "And generates tab-separated files formatted like ..."
	echo
	echo "output_index	spec_index	observed_MW	charge_state	scan_number	rank	calculated_MW	delta_mass	score	probability	peptide	protein	pept_position	mod1	mod2	mod3"
	echo "1	14	2869.4179	6	0	1	2870.6817	-1.2638	22	0.0004	R.A+5RVILGSPRPRVIVSSPWPAVVVASPR.P	sp|P0DI83	72~98	5	NA	NA"
	echo "1	14	2869.4179	6	0	2	2870.4205	-1.0026	18	0	R.NQCPAKPDGGGAPNGVRNGLAAELGPASPR.R	sp|Q9BX95	74~103	NA	NA	NA"
	echo "1	14	2869.4179	6	0	3	2868.3922	1.0257	15	0	R.QEPALRGSPGPLTPHPCNELGP+94PASPR.T	sp|Q8TAY7	39~65	94	NA	NA"
	echo "1	14	2869.4179	6	0	4	2869.4381	-0.0202	14	0	K.KKGPPPIR-5SCDDFSHMGTLPHSKSPR.Q	sp|O75815	55~80	-5	NA	NA"
	echo "1	14	2869.4179	6	0	5	2870.4587	-1.0408	13	0	R.PKAPASP-37RRPQTPTPSEQDADPGPASPR.D	sp|Q8NBB4	5~32	-37	NA	NA"
	echo "2	19	3094.4649	5	0	1	3094.465	-0.0001	82	0.9997	R.KTNDKDEKKEDGKQAENDSSNDDKTKK.S	sp|Q9BXP5	301~327	NA	NA	NA"
	echo "2	19	3094.4649	5	0	2	3094.3701	0.0949	75	0.9971	R.KTNDKDEKKEDGKQAENDSSNDDKT+128K.K	sp|Q9BXP5	301~326	128	NA	NA"
	echo "2	19	3094.4649	5	0	3	3094.3701	0.0949	71	0.9929	K.T+128NDKDEKKEDGKQAENDSSNDDKTKK.S	sp|Q9BXP5	302~327	128	NA	NA"
	echo "2	19	3094.4649	5	0	4	3094.5495	-0.0845	46	0.0083	K.VLL-33AWSGGPSSS+55SMVWQVLEGLSQDSAKR-2.L	sp|Q2VPK5	75~103	-33	55	-2"
	echo "2	19	3094.4649	5	0	5	3093.4716	0.9933	45	0.0137	K.KTL-41SQMSLSSDNSHATQNISP+78KKDD+50FK.N	sp|Q6ZP01	731~757	-41	78	50"

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





while [ $# -ne 0 ] ; do


#sed 's/\r/\n/g' $1 | awk '
perl -pe 's/\r/\n/g' $1 | awk '
function clear_variables(){
	output_index+=1;
	rank=0;
	spec_index="";
	observed_MW="";
	charge_state="";
	scan_number="";
}
BEGIN {
	FS="\t";
	OFS="\t";
	output_index=0;
	clear_variables();
	print "output_index	spec_index	observed_MW	charge_state	scan_number	rank	calculated_MW	delta_mass	score	probability	peptide	protein	pept_position	mod1	mod2	mod3"
}
( /^[[:blank:]]*$/ ) {
	clear_variables();
	next;
}
( $6 == "" && $7 == "" ){
	spec_index=$2;
	observed_MW=$3;
	charge_state=$4;
	scan_number=$5;
	next;
}
{
	rank+=1;
	calculated_MW=$1
	delta_mass=$2
	score=$3
	probability=$4
	peptide=$5
	protein=$6
	pept_position=$7
#	split($5,a,/[[:alpha:]+.]*/)
	split($5,a,/[A-Z+.]*/)
	mod1=(length(a) > 2) ? a[2] : "NA"
	mod2=(length(a) > 3) ? a[3] : "NA"
	mod3=(length(a) > 4) ? a[4] : "NA"
	print output_index,spec_index,observed_MW,charge_state,scan_number,rank,calculated_MW,delta_mass,score,probability,peptide,protein,pept_position,mod1,mod2,mod3
}
' > $1.OUTPUT.tsv

	shift

done
