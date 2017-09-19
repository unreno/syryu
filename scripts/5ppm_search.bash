#!/usr/bin/env bash

script=`basename $0`

#	Defaults:
verbose=false

function usage(){
	echo
	echo "Does stuff."
	echo
	echo "Usage:"
	echo
	echo "$script <OPTIONS> something_else"
	echo
	echo "Options:"
	echo "	--verbose ........... NO VALUE, just boolean flag"
	echo
	echo "Default option values:"
	echo "	--verbose .......... ${verbose}"
	echo
	echo "Examples:"
	echo "	$script"
	echo
	echo
	exit 1
}

while [ $# -ne 0 ] ; do
	case $1 in
#		-m|--m*)
#			shift; min=$1; shift;;
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


while [ $# -ne 0 ] ; do
	echo $1
	input=$1
	base=${1%.*}		#	drop the extension
	echo $base

#	
#	<RUN 2: 5ppm search>
#	Please use the same 3 mgf files in the RyuLab storage (PTM_blind).
#	 
#	1
#	Run MSGFPlus.jar with 3500M
#	Use Mods.txt for this search.
#	 


	java -Xmx3500M -jar /Users/jakewendt/Downloads/MSGFPlus/MSGFPlus.jar -s $input -o ${base}.mzid -d uniprot_reviewed_April2016.fasta -t 5ppm -ti 0,1 -inst 1 -protocol 1 -mod Mods.txt -n 10 -addFeatures 1
#	java -Xmx3500M -jar $(cygpath -w /cygdrive/c/ryulab/MSGFPlus/MSGFPlus.jar) -s $input -o ${base}.mzid -d uniprot_reviewed_April2016.fasta -t 5ppm -ti 0,1 -inst 1 -protocol 1 -mod Mods.txt -n 10 -addFeatures 1

#	2
#	Convert out file to tsv files (INPUTFILENAME is the same as OUTPUTFILENAME in step 2. >
#	 
	java -cp /Users/jakewendt/Downloads/MSGFPlus/MSGFPlus.jar edu.ucsd.msjava.ui.MzIDToTsv -i ${base}.mzid -showFormula 1
#	java -cp $(cygpath -w /cygdrive/c/ryulab/MSGFPlus/MSGFPlus.jar) edu.ucsd.msjava.ui.MzIDToTsv -i ${base}.mzid -showFormula 1
#	 

	shift
done
