#!/usr/bin/env bash

script=`basename $0`

function usage(){
	echo
	echo "Reads in given fasta file and counts amino acids in each sequence."
	echo "Output goes to filename with extension removed and '.amino_acid_counts.tsv' appended."
	echo
	echo "Usage:"
	echo
#	echo "$script <OPTIONS> fasta_files(s)"
	echo "$script fasta_files(s)"
	echo
#	echo "Options:"
#	echo "	--min INTEGER ......... Minimum value of shift range"
#	echo "	--verbose .................. NO VALUE, just boolean flag"
	echo
	echo "Default option values:"
#	echo "	--min ........ ${min}"
#	echo "	--verbose .......... ${verbose}"
	echo
	echo "Examples:"
	echo "	$script protein.fasta"
#	echo "	$script --min 75 protein.fasta"
	echo
	echo
	exit 1
}

while [ $# -ne 0 ] ; do
	case $1 in
#		-m|--m*)
#			shift; min=$1; shift;;
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
#	echo $1

	base=${1%.*}		#	drop the extension

	awk -v base=$base '
	function analyze_and_print() {
		gsub(/"/,"\"\"",description)
		printf( "\"%s\"\t",description ) >> out
		tmp=sequence
		for( i in acids )
			printf( "%s\t", gsub(acids[i],"",tmp) ) >> out
		print length(sequence) >> out
	}
	( !/^>/ ){
		gsub(/\r/,"");	#	remove windows carriage returns!
		sequence=sequence$0
	}
	( /^>/ ){
		if( description ) analyze_and_print()
		gsub(/^>/,"")
		description=$0
		sequence=""
	}
	BEGIN {
		out=base".amino_acid_counts.tsv"
		split("ABCDEFGHIJKLMNOPQRSTUVWXYZ",acids,"")
		printf("sequence\t") > out
		for( i in acids ) printf("%s\t",acids[i]) >> out
		printf("length\n") >> out
	}
	END{
		analyze_and_print()
	}' $1

	shift
done

