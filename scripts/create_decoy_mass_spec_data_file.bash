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
	echo "	--min_shift ........ Minimum value of shift range"
	echo "	--max_shift ........ Maximum value of shift range"
	echo "	--exclusion_min .... Do not shift values ABOVE this"
	echo "	--exclusion_max .... Do not shift values BELOW this"
	echo "	--number_of_peaks .. NOT IMPLEMENTED YET"
	echo "	--percent_of_peaks . PERCENT CHANCE PEAK WILL BE SHIFTED (not sure about this)"
	echo
	echo "Default option values:"
	echo "	--min_shift ........ ${min_shift}"
	echo "	--max_shift ........ ${max_shift}"
	echo "	--exclusion_min .... ${exclusion_min}"
	echo "	--exclusion_max .... ${exclusion_max}"
	echo "	--number_of_peaks .. ${number_of_peaks} NOT YET"
	echo "	--percent_of_peaks . ${percent_of_peaks}"
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
#					-b|--b*)
#						shift; blast=$1; shift ;;
		--min_shift)
			shift; min_shift=$1; shift;;
		--max_shift)
			shift; max_shift=$1; shift;;
		--exclusion_min)
			shift; exclusion_min=$1; shift;;
		--exclusion_max)
			shift; exclusion_max=$1; shift;;
		--percent_of_peaks)
			shift; percent_of_peaks=$1; shift;;
#number_of_peaks=
#percent_of_peaks=30
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
'
BEGIN{
	outfile=sprintf("%s.DECOY.mgf",infile);
	printf "" > outfile;
	srand();
}
( /^PEPMASS/ ){
	split($1,s,"=");
	pepmass = s[2];
	printf "PEPMASS = %s\n", pepmass;
}
( /^CHARGE/ ){
	split($1,s,"=");
	charge = s[2];
	printf "CHARGE = %s\n", charge;
	max_mass = pepmass * charge;
	printf "MAX MASS = %s\n", max_mass;
}
( /^[[:alpha:]]/ ){ print >> outfile; }
( /^[[:digit:]]/ ){
	print $1
	print exclusion_min, exclusion_max
	if( ( $1 < exclusion_min ) || ( $1 > exclusion_max ) ){
		print "... is not in exclusion range ...";
		if( int( rand() * 100 ) <= percent_of_peaks ) {
			print "... has been randomly selected to shift ...";
			random_shift = min_shift + int(rand() * ( max_shift - min_shift + 1 ) );
			new_mass = $1 + random_shift;
			print "... to new mass " new_mass " ...";
			if( ( new_mass > 0 ) && ( new_mass < max_mass ) ){
				print "... that is between 0 and " max_mass " ...";
				$1 = new_mass;
			} else {
				print "BUT IS NOT BETWEEN 0 and " max_mass;
			}
		}
	} else {
		print "... BUT IS IN EXCLUSION RANGE ...";
	}
	print $1, $2 >> outfile
}' $1

	shift

done
