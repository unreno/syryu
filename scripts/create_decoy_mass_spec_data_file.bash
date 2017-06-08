#!/usr/bin/env bash

script=`basename $0`
script_dir=`dirname $0`

#	Defaults:
min_shift=-200
max_shift=200
min_shifted_mass=0
exclusion_min=0
exclusion_max=0
number_of_peaks=
percent_of_peaks=30

function usage(){
	echo
	echo "Outputs file of same name with .DECOY.mgf appended to root name."
	echo
	echo "reads in an MGF file"
	echo "takes the following options"
	echo "	* min and max shift determines the range of shifting"
	echo "	* exclusion min and max determines range of masses to be left unshifted"
	echo "	* percent_of_peaks determines the chance of a given mass being shifted"
	echo "	* min shifted mass (default 0)"
	echo "reads its PEPMASS and CHARGE then computes a MAX SHIFTED MASS"
	echo "reads each MASS INTENSITY pair"
	echo "	* if mass is outside exclusion range"
	echo "		* percent chance of shifting determined"
	echo "	* and if shifting"
	echo "		* random number in shift range (positive and negative) determined"
	echo "	* and if summation of input mass and shift > MIN SHIFTED MASS and below MAX SHIFTED MASS"
	echo "		* mass is shifted"
	echo "	* else"
	echo "		* mass is left as is"
	echo "outputs input MGF file with (possibly) modified mass values"
	echo
	echo "Usage:"
	echo
	echo "$script <OPTIONS> mgf_file(s)"
	echo
	echo "Options:"
	echo "	--min_shift INTEGER......... Minimum value of shift range"
	echo "	--max_shift INTEGER......... Maximum value of shift range"
	echo "	--min_shifted_mass INTEGER.. Minimum value for shifted mass"
	echo "	--exclusion_min INTEGER..... Do not shift values ABOVE this"
	echo "	--exclusion_max INTEGER..... Do not shift values BELOW this"
	echo "	--number_of_peaks .......... NOT IMPLEMENTED YET"
	echo "	--percent_of_peaks INTEGER.. PERCENT CHANCE PEAK WILL BE SHIFTED (not sure about this)"
	echo "	--verbose .................. NO VALUE, just boolean flag."
	echo "	          .................. Output file will be invalid MGF as contains a lot of other stuff"
	echo
	echo "Default option values:"
	echo "	--min_shift ......... ${min_shift}"
	echo "	--max_shift ......... ${max_shift}"
	echo "	--min_shifted_mass .. ${min_shifted_mass}"
	echo "	--exclusion_min ..... ${exclusion_min}"
	echo "	--exclusion_max ..... ${exclusion_max}"
	echo "	--number_of_peaks ... ${number_of_peaks} NOT YET"
	echo "	--percent_of_peaks .. ${percent_of_peaks}"
	echo "	--verbose ........... ${verbose}"
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

while [ $# -ne 0 ] ; do
	case $1 in
		--min_shift)
			shift; min_shift=$1; shift;;
		--max_shift)
			shift; max_shift=$1; shift;;
		--min_shifted_mass)
			shift; min_shifted_mass=$1; shift;;
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

#	Basically, this is TRUE AND DO ...
[ $# -eq 0 ] && usage

#number_of_peaks=
#percent_of_peaks=30

while [ $# -ne 0 ] ; do


#	outfile=sprintf("%s.DECOY.mgf",infile_base);
#	printf "" > outfile;
#( /^[[:alpha:]]/ ){ print >> outfile; }
#	print $1, $2 >> outfile

awk \
	-v min_shift=$min_shift \
	-v max_shift=$max_shift \
	-v min_shifted_mass=$min_shifted_mass \
	-v exclusion_min=$exclusion_min \
	-v exclusion_max=$exclusion_max \
	-v percent_of_peaks=$percent_of_peaks \
	-v verbose=$verbose \
'
BEGIN{ srand(); }
( /^PEPMASS/ ){
	split($1,s,"=");
	pepmass = s[2];
	if( verbose ) printf( "PEPMASS = %s\n", pepmass );
}
( /^CHARGE/ ){
	split($1,s,"=");
	charge = s[2];
	if( verbose ) printf( "CHARGE = %s\n", charge );
	max_shifted_mass = pepmass * charge;
	if( verbose ) printf( "MAX MASS = %s\n", max_shifted_mass );
}
( /^[[:alpha:]]/ ){ print }
( /^[[:digit:]]/ ){
	if( verbose ) print $1
	if( verbose ) print "Exclusion range " exclusion_min " - " exclusion_max
	if( ( $1 < exclusion_min ) || ( $1 > exclusion_max ) ){
		if( verbose ) print "... is not in exclusion range ...";
		num = int( rand() * 100 )
		if( num < percent_of_peaks ) {
			if( verbose ) print "... has been randomly ( " num " < " percent_of_peaks " ) selected to shift ...";
			if( verbose ) print "... between " min_shift " and " max_shift " ...";

			original_mass = $1;
			i=0

			#	try 3 times to shift, if not then just skip.
			while( original_mass == $1 && i < 3 ){
				if( verbose && i > 0 ) print "TRYING AGAIN";

				random_shift = min_shift + int(10 * rand() * ( max_shift - min_shift + 1 ) ) / 10;
				if( verbose ) print "... by " random_shift " ...";
				new_mass = $1 + random_shift;
				if( verbose ) print "... to new mass " new_mass " ...";
				if( ( new_mass > min_shifted_mass ) && ( new_mass < max_shifted_mass ) ){
					if( verbose ) print "... that is between " min_shifted_mass " and " max_shifted_mass " ...";
					$1 = new_mass;
				} else {
					if( verbose ) print "BUT IS NOT BETWEEN " min_shifted_mass " and " max_shifted_mass;
				}
				i++;
			}

		}
	} else {
		if( verbose ) print "... BUT IS IN EXCLUSION RANGE ...";
	}
	print $1, $2
}' $1 | awk '
BEGIN { split("",buffer) }
( /^[0-9 .]+$/ ){
	buffer[length(buffer)+1]=$0
}
( !/^[0-9 .]+$/ ){
	if( length(buffer) > 0 ){
		asort(buffer,sorted,"@val_num_asc")
		for(i=1;i<=length(sorted);i++)
			print sorted[i];
		delete(buffer);
	}
	print
}' > ${1%.*}.DECOY.mgf

	shift

done


#	Create a decoy mass spectrometry data file
#
#	Goal: to construct a null distribution of PTM matching scores
#	How to:
#	Shift masses of some peaks randomly in a given mgf file
#
#
#	Exclude some peaks in a certain mass range to be shifted
#	For example,
#		if the exclusion mass range is between 0 and 200,
#		then any masses between this range cannot be shifted.
#
#	Minimum mass shift and Maximum mass shift
#	For example,
#		if a minimum mass shift is 10 and a maximum mass shift is 200,
#		then a random number from (10, 200) will be chosen.
#		If 5.5 was chosen for a certain peak with mass 130.5,
#		then the shifted mass for that peak is 136.
#
#	Option to choose either number of peaks to be shifted or % of peaks to be shifted
#
#	Input file (mgf file)
#	Output file  (mgf file with another name)
#
#
#	Restriction
#
#	The minimum mass cannot be less than zero.
#	The maximum mass cannot be larger than Pepmass*Charge
#
#
#	Example
#
#	The file contains 2 spectra.
#	Spectra 1 (PepMass 1000 and charge 1)
#	The following is masses.
#	135.5
#	143.4 -> +5.5 =>148.9	(random chosen to shift)
#	200.5
#	300.4 -> -10 => 290.4	(randomly chosen to shift)
#	420.6
#	900.5
#	Spectra 2 (PepMass 500 and charge 2)
#	300.2 -> 50 => 350.2	(randomly chosen to shift)
#	403.50
#	950.5 -> +60 =>  (randomly chosen to shift however CANNOT DO THIS BECAUSE IT BECOME GREATED THAN 500*2)
#



#	EXAMPLE SPECTRUM
#	Only PEPMASS, CHARGE and MASS (first column) values are considered in this script


#	BEGIN IONS
#	TITLE=140521_EOC_MCis_NT_3.12.12.2 File:"140521_EOC_MCis_NT_3.raw", NativeID:"controllerType=0 controllerNumber=1 scan=12"
#	RTINSECONDS=3.41647884
#	PEPMASS=573.837194615683 4476.82373046875			[second value not relevant here]
#	CHARGE=2+
#	148.3734771 544.1920166016				[line in format of MASS INTENSITY]
#	149.0221207 11944.6396484375
#	150.0263882 509.1209106445
#	151.1748135 597.6806030273
#	158.5967596 557.734375
#	161.7232055 524.9778442383
#	178.6890678 590.6413574219
#	195.6584047 499.9047241211
#	205.0836878 1168.0891113281
#	293.1799579 564.6156616211
#	307.3158042 567.4732055664
#	845.83007 616.4040527344
#	926.6448323 681.1176757813
#	930.4296215 634.0212402344
#	END IONS




