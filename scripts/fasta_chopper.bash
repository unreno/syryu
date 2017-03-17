#!/usr/bin/env bash


script=`basename $0`

#	Defaults:
min=50

function usage(){
	echo
	echo "Splits reads in given fasta file into sequences at given min length."
	echo
	echo "Usage:"
	echo
	echo "$script <OPTIONS> fasta_files(s)"
	echo
	echo "Options:"
	echo "	--min INTEGER ......... Minimum value of shift range"
#	echo "	--verbose .................. NO VALUE, just boolean flag"
	echo
	echo "Default option values:"
	echo "	--min ........ ${min}"
#	echo "	--verbose .......... ${verbose}"
	echo
	echo "Examples:"
	echo "	$script protein.fasta"
	echo "	$script --min 75 protein.fasta"
	echo
	echo
	exit 1
}

while [ $# -ne 0 ] ; do
	case $1 in
		-m|--m*)
			shift; min=$1; shift;;
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

	awk -v min=$min '
	function split_and_print() {
		i=1
		name1=description
		name2=""
		if( ( space_position = index(description," ") ) > 0 ){
			name1=substr(description,0,space_position-1)
			name2=substr(description,space_position)
		}
		while( ( l = length(sequence) ) > 0 ){
			print name1 "_" i name2
			if( l < ( 2 * min ) ) {
				print sequence
				sequence=""
			} else {
				print substr(sequence, 0, min )
				i+=1
				sequence=substr(sequence,min+1)
			}
		}
	
	}
	( !/^>/ ){
		gsub(/\r/,"");	#	remove windows carriage returns!
		sequence=sequence$0
	}
	( /^>/ ){
		if( description ) split_and_print()
		description=$0
		sequence=""
	}
	END{
		split_and_print()
	}' $1

	shift
done

