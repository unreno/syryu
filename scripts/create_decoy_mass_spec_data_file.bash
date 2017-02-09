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
  echo
  echo "Usage:"
  echo
  echo "$script <OPTIONS> mgf_file(s)"
  echo
#			  echo "COMMANDS: blastn, tblastx"
	echo "Options:"
	echo "	--min_shift ........ Minimum value of shift range"
	echo "	--max_shift ........ Maximum value of shift range"
	echo "	--exclusion_min .... ${exclusion_min}"
	echo "	--exclusion_max .... ${exclusion_max}"
	echo "	--number_of_peaks .. ${number_of_peaks}"
	echo "	--percent_of_peaks . ${percent_of_peaks}"
	echo
	echo "Default option values:"
	echo "	--min_shift ........ ${min_shift}"
	echo "	--max_shift ........ ${max_shift}"
	echo "	--exclusion_min .... ${exclusion_min}"
	echo "	--exclusion_max .... ${exclusion_max}"
	echo "	--number_of_peaks .. ${number_of_peaks}"
	echo "	--percent_of_peaks . ${percent_of_peaks}"
	echo
	echo "Notes:"
	echo "	number and percent peaks are mutually exclusive. Setting one deletes the other."
	echo "	options are specified WITHOUT EQUALS SIGN. (ie. --min_shift -50 )"
  echo
  echo "Example:"
#			  echo "$script dna/output/fallon_SFPB001A_filtered_20130722/trinity_input_single.fasta.blastn.txt"
  echo
	exit 1
}
#	Basically, this is TRUE AND DO ...
[ $# -eq 0 ] && usage

while [ $# -ne 0 ] ; do
	case $1 in
#					-b|--b*)
#						shift; blast=$1; shift ;;
#min_shift=-200
#max_shift=200
#exclusion_min=0
#exclusion_max=0
#number_of_peaks=
#percent_of_peaks=30
		-*)
			echo ; echo "Unexpected args from: ${*}"; usage ;;
		*)
			break;;
	esac
done


exit


	
#	
#	Create a decoy mass spectrometry data file
#	Goal: to construct a null distribution of PTM matching scores
#	How to:
#	Shift masses of some peaks randomly in a given mgf file
#	
#	
#	Parameters
#	Exclude some peaks in a certain mass range to be shifted
#	For example, if the exclusion mass range is between 0 and 200, then any masses between this range cannot be shifted.
#	Minimum mass shift and Maximum mass shift
#	For example, if a minimum mass shift is 10 and a maximum mass shift is 200, then a random number from (10, 200) will be chosen. If 5.5 was chosen for a certain peak with mass 130.5, then the shifted mass for that peak is 136. 
#	Option to choose either number of peaks to be shifted or % of peaks to be shifted
#	Input file (mgf file)
#	Output file  (mgf file with another name)
#	
#	Restriction
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
#	143.4 -> +5.5 =>148.9		#	( 143.5 + random number( 5.5 ) )
#	200.5
#	300.4 -> -10 => 290.4	#	(300.4 + random number (-10) )
#	420.6
#	900.5
#	Spectra 2 (PepMass 500 and charge 2)
#	300.2 -> 50 => 350.2	#	 300.2 + 50
#	403.50
#	950.5 -> +60 =>  (CANNOT DO THIS BECAUSE IT BECOME GREATER THAN 500*2)
#	
#	
#	
#	BEGIN IONS
#	TITLE=140521_EOC_MCis_NT_3.12.12.2 File:"140521_EOC_MCis_NT_3.raw", NativeID:"controllerType=0 controllerNumber=1 scan=12"
#	RTINSECONDS=3.41647884
#	PEPMASS=573.837194615683 4476.82373046875
#	CHARGE=2+
#	148.3734771 544.1920166016
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
#	
#	
#	PEPMASS=INTENSITY .... irrelevant
#	CHARGE=CHARGE ... self explanatory
#	
#	First number is MASS, Second is INTENSITY
#	
#	
#	
#	
#	


##	blast_check.gawk
#			BEGIN{
#				if( head == "" ) head="BLASTN"
#				if( tail == "" ) tail="Gap Penalties"
#				head_count=0
#				tail_count=0
#			}
#			(NR%100000 == 0){
#				print "Read",FNR,"lines from",FILENAME >> "blast_check.log"
#			}
#			( $0 ~ head ){ head_count++ }
#			( $0 ~ tail ){ tail_count++ }
#			(/Query= /){ query++ }
#			(/Effective search space used:/){ effective++ }
#			(/no longer exists in database/){ nolonger++ }
#			(/[^[:print:]]/){ 	#	too many (includes \anything?)
#				print "Non-printable character(s) found on line number :",NR,":"
#				print $0
#				nonprint++ 
#			}
#			#	control chars found on EVERY line if this were awk instead of gawk
#			(/[[:cntrl:]]/){	#	REQUIRES double brackets as they are for different reasons
#				print "Control character(s) found on line number :",NR,":"
#				print $0
#				control++ 
#			}
#			END{
#				print "First line count :",head_count,":"
#				print "Last line count :",tail_count,":"
#				if( head_count != tail_count ){ print " * First and Lines lines are out of sync" }
#				print "Query=  line count :",query,":"
#				print "Effective search space used: line count :",effective,":"
#				if( query != effective ){ print " * Query= and Effective search lines are out of sync" }
#				print "no longer exists in database line count :",nolonger,":"
#				if( nolonger > 0 ){ print " * TOO MANY no longer exists in database lines" }
#				print "nonprint character line count :",nonprint,":"
#				if( nonprint > 0 ){ print " * TOO MANY nonprintable characters" }
#				print "control character line count :",control,":"
#				if( control > 0 ){ print " * TOO MANY control characters" }
#				print "---"
#				close("blast_check.log")
#			}

