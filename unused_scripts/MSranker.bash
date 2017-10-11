#!/usr/bin/env bash


script=`basename $0`
script_dir=`dirname $0`

#	Defaults:


function usage(){
	echo
	echo $script
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
	echo "B.   Using “theospec” codes (C codes), calculate theoretical m/z value for the sequence. “theospec” has a slightly different peptide sequence format than MODa output."
	echo "You need to either modify our peptide sequence in MODa or modify theospec codes."
	echo "What theospec can do:"
	echo "For the first sequence, I need a command “./theospec -z5 R.KTNDKDEKKEDGKQAENDSSNDDKTKK.S”"
	echo "For the second sequence, I need a command “./theospec -z5 -c*=128 R.KTNDKDEKKEDGKQAENDSSNDDKT*K.K”"
	echo "We may have more than one modification, so if we have “R.KTND-12KDEKKEDGKQAENDSSNDDKT+128K.K”,"
	echo "The command is like this: ./theospec -z1 -c*=128 -c#=-12 R.KTND#KDEKKEDGKQAENDSSNDDKT*K.K"
	echo
	echo "For now, printing output in the output file will do. But in the future, you will pass these output values to another function that you will write (#4 step)."
	echo "I attached out_blind_140521….txt.OUTPUT."
	echo "Theospec is in https://sourceforge.net/projects/protms/files/theospec/3.00/"
	echo
	echo
	echo
	echo "The step #3, #4, and #5 are very complicated. It is the following:"
	echo "a.       For theoretical spectra (i.e. 1535.1.theospec.txt), locate the corresponding scan in mgf file (i.e. 140521_EOC_MCis_T2_3.SCANS.mgf). In Azure storage, there are mgf files. Please use 140521_EOC_MCis_T2_3.SCANS.mgf."
	echo "For example, for 1535.1.theospec.txt, the scan number is 1535. (Note that we have multiple results for the same scan. 1535.2, 1535.3, and so on.)."
	echo "In mgf, the spectrum with SCANS=1535 is the one we want for 1535.1.theospec.txt."
	echo
	echo "b.        In the located mgf spectrum, there are m/z values (1st column) and intensities (2nd column)."
	echo "Goal: Find the matched peaks by their m/z values with given tolerance (let’s call it “tol”. Let use a default tol=0.5.)."
	echo "For each m/z value in corresponding spectrum in mgf (denoted as observed m/z),"
	echo "look for whether there exists (a) m/z value(s) in theoretical spectrum (denoted as theoretical m/z) within the specified +/- tol. For example, if 200 m/z value is in mgf file, any m/z values between [200-tol, 200+tol] are potential matches."
	echo
	echo "The first column is m/z value in theospec,"
	echo "From theospec,"
	echo "106.5733148696(65);+2;b2;0;(-H2O)"
	echo "109.7269497461(68);+3;b3;0;(-NH3)"
	echo "111.0577632881(70);+4;b4;0;(-H2O)"
	echo "112.0756903831(62);+1;b1;0;(-NH3)"
	echo "115.4024661132(75);+3;b3;0;()"
	echo "115.5604044596(76);+4;b4;0;()"
	echo "115.5785972128(76);+2;b2;0;()"
	echo "…"
	echo
	echo "Then, in this theospec, 106.5733148696, 109.7269497461, … are m/z values."
	echo
	echo "We need to have one-to-one matching between observed and theoretical m/z values."
	echo "If there are more than one matched theoretical m/z per one observed m/z, then pick the closest theoretical m/z value to the observed one. This observed peak is matched."
	echo "If there are more than one matched observed m/z per theoretical m/z, then pick the m/z values with the largest intensities. This theoretical peak is matched."
	echo
	echo
	echo "c.       For each pair of theospec file & matched mgf spectrum, please record the following:"
	echo "-The total number of observed m/z values"
	echo "-The total number of theoretical m/z values"
	echo "-The largest intensity in spectrum whether it is matched or not."
	echo "-# of observed b-ion m/z values matched (b-ion is starting with b in theoretical spec; for example, b2, b3, b4, ….) Let’s call it as “Nb”"
	echo "-# of observed y-ion m/z values matched (y-ion is starting with y in theoretical spec; for example, y2, y5,…) Let’s call it as “Ny”"
	echo "-The list of matched observed b-ion intensities in the spectrum (in mgf files). – separated by semicolon"
	echo "-The list of matched observed b-ion intensities divided by the largest intensity in spectrum, and times 100. Let’s call it “Intensity_b_i” for each ith observed peak. – separated by semicolon"
	echo "- The list of matched observed y-ion intensities (in mgf files) – separated by semicolon"
	echo "-The list of matched observed y-ion intensities divided by the largest intensity in spectrum, and times 100. Let’s call it “Intensity_y_i” for each ith observed peak. – separated by semicolon"
	echo "-m/z errors between observed and theoretical m/z values for b-ion matches – separated by semicolon"
	echo "-m/z errors between observed and theoretical m/z values for y-ion matches – separated by semicolon"
	echo "                Intermediate file can ends with “ion_statistics”"
	echo
	echo "d.       We will try various way to calculated the matching scores between observed and theoretical m/z values using ion_statistics. However, this will be the first scoring (let’s call it hyperscore)."
	echo "For this formula, we will use only underlined values in step c."
	echo "The formula is the following:"
	echo "LOG (Nb! * Ny! * (Sum of intensity_b_i from i=1 to Nb) * (Sum of intensity_y_i from i=1 to Ny))."
	echo "Here is log is logarithm with base e."
	echo "! is a factorial."
	echo
	echo "For example,"
	echo "Assume that we have 3 matched b-ions. And their intensity_b_i values were 2, 50, 4."
	echo "Assume that we have 4 matched y-ions. And their intensity_y_i values were 1,4,100,30."
	echo "Then, log( 3! * 4! * (2+50+4) * (1+ 4+100+300)). The answer is 14.99905."
	echo
	echo "e.       Now, we go back to our original file, out_blind_140521….txt.OUTPUT and attach this hyperscore in each corresponding line (matched with scan number and rank)."
	echo
	echo
	echo "Usage:"
	echo
	echo "$script <OPTIONS> input_file mgf_spectrum_file"	#(s)"
	echo
	echo "Will create a folder containing all of the theospec call output."
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
[ $# -ne 2 ] && usage

date=$(date "+%Y%m%d%H%M%S")


date="TESTING"



#
#	Multple files are passed for each processing so DO NOT loop over them.
#
#while [ $# -ne 0 ] ; do

theodir="${1}_theospec_${date}/"
obsdir="${1}_obsspec_${date}/"

if [ ! -d $theodir ] ; then

	mkdir -p "$theodir"

	#	gawk -F"\t" '
	#	MUST use gawk to use join
	perl -pe 's/\r//g' $1 | gawk -v dir=$theodir '
	@include "join"
	BEGIN{ 
		FS="\t"
		OFS="\t"
	}
	( NR == 1 ){ next }
	{
		split($11,parts,/\./)		#	-.QWLPKD-123EESFLQFKYQALQVP+117.P	(theoretically could include decimals)
		a=parts[1]							#	-
		b=join(parts,2,length(parts)-1,".")	#	QWLPKD-123EESFLQFKYQALQVP+117
		c=parts[length(parts)]	#	P

		#	parse middle extracting numbers and replacing them with codes passed to theospec
		split(b,numbers,/[A-Z]+/)	#	["",-123,+117]
		split(b,letters,/[0-9.+-]+/)	#	[	"QWLPKD","EESFLQFKYQALQVP",""]

		outpeptide=""
		for(i=1;i<=length(letters);i++){
			outpeptide=sprintf("%s%s",outpeptide,letters[i])
			if( length(numbers[i+1]) > 0 )
				outpeptide=sprintf("%s%d",outpeptide,i)
		}
		#	"QWLPKD1EESFLQFKYQALQVP2"

		z=( $4 > 5 ) ? 5 : $4

		theospec_command="theospec -z"z
		for(i=2;i<length(numbers);i++){
			theospec_command=sprintf("%s -c%d=%s",theospec_command,i-1,numbers[i])
		}
		#	-c1=-123 -c2=+117

		theospec_command=sprintf("%s %s.%s.%s > %s%s.%s.theospec.txt",theospec_command,a,outpeptide,c,dir,$5,$6)
		print theospec_command
	}' | tee $theodir/theospec_command_list.txt | bash

fi


# > $1.peptide_and_charge.tsv
#' > $1.MSranker.tsv









#	For the final output file, appending MSranker.tsv will be good. For intermediate output files for step 1 and 2, filename_theoSpec.tsv will be good.


#	#	remove any empty error files
#	for f in $(find $dir -name *errors) ; do
#		if [ ! -s $f ] ; then
#			rm $f
#		fi
#	done





if [ ! -d $obsdir ] ; then

	mkdir -p "$obsdir"

	awk -v scan=$scan -v dir=$obsdir '
		( /^BEGIN/ ){
			buffer=$0; next;
		}
		( match( $0, "SCANS=" ) ){
			split($0,a,"=")
			scan=a[2];
		}
		{
			buffer=buffer "\n" $0
		}
		( /^END/ ){
			print buffer > dir"/"scan".observed.txt"
		}
	' $2


fi



for f in $(find $theodir -name *theospec.txt) ; do
	#echo $f
	base=$(basename $f)
	scan=${base%%.*}
	#echo $scan

	obsfile="$obsdir/$scan.observed.txt"

	if [ ! -f $obsfile ] ; then
		echo "Missing $obsfile"
	else

		#	loop over observed file m/z / intensities

		awk '( /^[[:digit:]]/ )' $obsfile | while read -r $mzi ; do

			#	find matching theospec m/z's ???
			#	Goal: Find the matched peaks by their m/z values with given tolerance (let’s call it “tol”. Let use a default tol=0.5.).
			#	For each m/z value in corresponding spectrum in mgf (denoted as observed m/z),
			#	look for whether there exists (a) m/z value(s) in theoretical spectrum (denoted as theoretical m/z) within the specified +/- tol. For example, if 200 m/z value is in mgf file, any m/z values between [200-tol, 200+tol] are potential matches.


			awk -v obsmzi=$mzi 'BEGIN{
				split(obsmzi,obs)
			}{
				print obsmzi 
			}' $f

		done

	fi

done

#' > $1.MSranker.tsv
