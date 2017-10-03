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


MSGFPlus=MSGFPlus.jar

echo $OSTYPE 
case "$OSTYPE" in
#	solaris*) echo "SOLARIS" ;;
	darwin*)  MSGFPlus=/Users/jakewendt/Downloads/MSGFPlus/MSGFPlus.jar ;;
#	linux*)   echo "LINUX" ;;
#	bsd*)     echo "BSD" ;;
#	msys*)    echo "WINDOWS" ;;
	cygwin*)  MSGFPlus=$(cygpath -w /cygdrive/c/ryulab/MSGFPlus/MSGFPlus.jar) ;;
#	*)        echo "unknown: $OSTYPE" ;;
esac


while [ $# -ne 0 ] ; do
	echo $1
	input=$1
	base=${1%.*}		#	drop the extension
	echo $base

	cmd="java -Xmx3500M -jar $MSGFPlus -s $input -o ${base}.5ppm.mzid -d uniprot_reviewed_April2016.fasta -t 5ppm -ti 0,1 -inst 1 -protocol 1 -mod Mods.txt -n 10 -addFeatures 1"
	echo $cmd
	$cmd

	cmd="java -cp $MSGFPlus edu.ucsd.msjava.ui.MzIDToTsv -i ${base}.5ppm.mzid -showFormula 1"
	echo $cmd
	$cmd

	shift
done
