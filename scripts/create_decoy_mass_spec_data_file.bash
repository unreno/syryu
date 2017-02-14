#!/usr/bin/env bash

script=`basename $0`

#	Defaults:
min_shift=-200
max_shift=200
exclusion_min=0
exclusion_max=0
number_of_peaks=
percent_of_peaks=30

function usage(){
  echo
#			  echo "checks the contents of the given blast output file"
	echo "Outputs file file of same name with .DECOY.mgf appended to root name."
  echo
  echo "Usage:"
  echo
  echo "$script <OPTIONS> mgf_file(s)"
  echo
	echo "Options:"
	echo "	--min_shift INTEGER......... Minimum value of shift range"
	echo "	--max_shift INTEGER......... Maximum value of shift range"
	echo "	--exclusion_min INTEGER..... Do not shift values ABOVE this"
	echo "	--exclusion_max INTEGER..... Do not shift values BELOW this"
	echo "	--number_of_peaks .......... NOT IMPLEMENTED YET"
	echo "	--percent_of_peaks INTEGER.. PERCENT CHANCE PEAK WILL BE SHIFTED (not sure about this)"
	echo "	--verbose .................. NO VALUE, just boolean flag"
	echo
	echo "Default option values:"
	echo "	--min_shift ........ ${min_shift}"
	echo "	--max_shift ........ ${max_shift}"
	echo "	--exclusion_min .... ${exclusion_min}"
	echo "	--exclusion_max .... ${exclusion_max}"
	echo "	--number_of_peaks .. ${number_of_peaks} NOT YET"
	echo "	--percent_of_peaks . ${percent_of_peaks}"
	echo "	--verbose .......... ${verbose}"
	echo
	echo "Notes:"
	echo "	number and percent peaks are mutually exclusive. Setting one deletes the other."
	echo "	options are specified WITHOUT EQUALS SIGN. (ie. --min_shift -50 )"
  echo
  echo "Examples:"
	echo "	$script 140521_EOC_MCis_NT_3-head-20.mgf"
	echo "	$script --exclusion_max 900 --percent_of_peaks 100 140521_EOC_MCis_NT_3.mgf"
	echo
  echo
	exit 1
}
#	Basically, this is TRUE AND DO ...
[ $# -eq 0 ] && usage

while [ $# -ne 0 ] ; do
	case $1 in
		--min_shift)
			shift; min_shift=$1; shift;;
		--max_shift)
			shift; max_shift=$1; shift;;
		--exclusion_min)
			shift; exclusion_min=$1; shift;;
		--exclusion_max)
			shift; exclusion_max=$1; shift;;
		--percent*)
			shift; percent_of_peaks=$1; shift;;
		-v|--v*)
			verbose=true; shift;;
		-*)
			echo ; echo "Unexpected args from: ${*}"; usage ;;
		*)
			break;;
	esac
done

#number_of_peaks=
#percent_of_peaks=30

while [ $# -ne 0 ] ; do


awk \
	-v infile=$1 \
	-v min_shift=$min_shift \
	-v max_shift=$max_shift \
	-v exclusion_min=$exclusion_min \
	-v exclusion_max=$exclusion_max \
	-v percent_of_peaks=$percent_of_peaks \
	-v verbose=$verbose \
'
BEGIN{
	outfile=sprintf("%s.DECOY.mgf",infile);
	printf "" > outfile;
	srand();
}
( /^PEPMASS/ ){
	split($1,s,"=");
	pepmass = s[2];
	if( verbose ) printf( "PEPMASS = %s\n", pepmass );
}
( /^CHARGE/ ){
	split($1,s,"=");
	charge = s[2];
	if( verbose ) printf( "CHARGE = %s\n", charge );
	max_mass = pepmass * charge;
	if( verbose ) printf( "MAX MASS = %s\n", max_mass );
}
( /^[[:alpha:]]/ ){ print >> outfile; }
( /^[[:digit:]]/ ){
	if( verbose ) print $1
	if( verbose ) print "Exclusion range " exclusion_min " - " exclusion_max
	if( ( $1 < exclusion_min ) || ( $1 > exclusion_max ) ){
		if( verbose ) print "... is not in exclusion range ...";
		num = int( rand() * 100 )
		if( num < percent_of_peaks ) {
			if( verbose ) print "... has been randomly ( " num " < " percent_of_peaks " ) selected to shift ...";
			if( verbose ) print "... between " min_shift " and " max_shift " ...";
			random_shift = min_shift + int(10 * rand() * ( max_shift - min_shift + 1 ) ) / 10;
			if( verbose ) print "... by " random_shift " ...";
			new_mass = $1 + random_shift;
			if( verbose ) print "... to new mass " new_mass " ...";
			if( ( new_mass > 0 ) && ( new_mass < max_mass ) ){
				if( verbose ) print "... that is between 0 and " max_mass " ...";
				$1 = new_mass;
			} else {
				if( verbose ) print "BUT IS NOT BETWEEN 0 and " max_mass;
			}
		}
	} else {
		if( verbose ) print "... BUT IS IN EXCLUSION RANGE ...";
	}
	print $1, $2 >> outfile
}' $1

	shift

done