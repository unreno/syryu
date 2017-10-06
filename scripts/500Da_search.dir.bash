#!/usr/bin/env bash

script=`basename $0`

#	Defaults:
verbose=false

function usage(){
	echo
	echo "Runs MSGFPlus on mgf files in given directory."
	echo
	echo "Usage:"
	echo
	echo "$script <OPTIONS> <directory containing mgf files>"
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

#MSGFPlus=/Users/jakewendt/Downloads/MSGFPlus/MSGFPlus.jar
#MSGFPlus=$(cygpath -w /cygdrive/c/ryulab/MSGFPlus/MSGFPlus.jar)


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

#	split_mgf_dir=$(mgf_splitter.bash --max 10 $1)

	for file in $( find $1 -name \*mgf	); do
		
		base=${file%.*}		#	drop the extension

		cmd="java -Xmx16G -jar $MSGFPlus -s $file -o $base.500Da.mzid -d uniprot_reviewed_April2016.fasta -t 500Da -ti 0,1 -inst 1 -protocol 1 -mod Mods2.txt -n 10 -addFeatures 1"
		echo $cmd
		$cmd
		if [ $? -ne 0 ] ; then
			echo $cmd >> failed_commands
		fi

#	check return status?
#	check presence of expected output file?
#	Do what if fails?


		cmd="java -cp $MSGFPlus edu.ucsd.msjava.ui.MzIDToTsv -i $base.500Da.mzid -showFormula 1"
		echo $cmd
		$cmd
		if [ $? -ne 0 ] ; then
			echo $cmd >> failed_commands
		fi



#	check return status?
#	check presence of expected output file?
#	Do what if fails?





	done






	base=${1%.*}		#	drop the extension
#	echo $base

	#	Merge resulting tsv files

#	first_tsv=$(ls -1 $split_mgf_dir/*tsv | head -1 )
#	head -1 $first_tsv > $base.500Da.tsv
#	tail --quiet --lines=+2 $split_mgf_dir/*tsv >> $base.500Da.tsv


#	/bin/rm -rf $split_mgf_dir



	shift
done
