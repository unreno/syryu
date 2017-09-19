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

	base=${1%.*}		#	drop the extension
	echo $base


#	
#	<RUN 1: 500Da search>
#	I uploaded 3 mgf files in the RyuLab storage (PTM_blind).
#	 
#	1
#	Break mgf files into smaller files (5 scans per file)
#	 
#	ex) The beginning of scan is "BEGIN IONS" and the end of scan is "END IONS".
#	 
#	 
#	2
#	Run MSGFPlus.jar with 10000M
#	Use Mods2.txt instead of Mods.txt for this search.
#	 
	java -Xmx10000M -jar /Users/jakewendt/Downloads/MSGFPlus/MSGFPlus.jar -s $1 -o $base.mzid -d uniprot_reviewed_April2016.fasta -t 500Da -ti 0,1 -inst 1 -protocol 1 -mod Mods2.txt -n 10 -addFeatures 1
#	java -Xmx10000M -jar $(cygpath -w /cygdrive/c/ryulab/MSGFPlus/MSGFPlus.jar) -s <INPUTFILENAME>.mgf -o <OUTPUTFILENAME>.mzid -d uniprot_reviewed_April2016.fasta -t 500Da -ti 0,1 -inst 1 -protocol 1 -mod Mods2.txt -n 10 -addFeatures 1
#	 
#	 
#	3
#	Convert output file to tsv files (INPUTFILENAME is the same as OUTPUTFILENAME in step 2. >
#	 
#	java -cp $(cygpath -w /cygdrive/c/ryulab/MSGFPlus/MSGFPlus.jar) edu.ucsd.msjava.ui.MzIDToTsv  -i  <INPUTFILENAME>.mzid -showFormula 1
#	 
#	 
#	4
#	Merge resulting tsv files
#	 

	shift
done
