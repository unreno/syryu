#!/usr/bin/env bash


script=`basename $0`

#	Defaults:
max=5

function usage(){
	echo
	echo "Splits mgf file limiting the number of scans."
	echo "Output goes to filename with number inserted between base and extension" 
	echo
	echo "Usage:"
	echo
	echo "$script <OPTIONS> mgf_files(s)"
	echo
	echo "Options:"
	echo "	--max INTEGER ......... Maximum number of scans per file"
#	echo "	--verbose .................. NO VALUE, just boolean flag"
	echo
	echo "Default option values:"
	echo "	--max ........ ${max}"
#	echo "	--verbose .......... ${verbose}"
	echo
	echo "Examples:"
	echo "	$script file.mgf"
	echo "	$script --max 10 file.mgf"
	echo
	echo
	exit 1
}

while [ $# -ne 0 ] ; do
	case $1 in
		-m|--m*)
			shift; max=$1; shift;;
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
	ext=${1##*.}		#	drop the base

#	echo $base
#	echo $ext

	basename=$(basename $base)

	scan_count=$(grep -c "^BEGIN IONS" $1)
#	echo $scan_count

	date=$(date "+%Y%m%d%H%M%S")
#	echo $date

	dir="${base}_${date}"
	mkdir -p "$dir"

	awk -v max=$max -v dir=$dir -v basename=$basename -v ext=$ext -v scan_count=$scan_count '
	BEGIN {
		i=0
		scans=0
	}
	( /^BEGIN IONS/ ){
		scans++
	}
	( /^BEGIN IONS/ && ( scans >= max ) ){
		scans=0
		i++
	}
	{
		print >> dir"/"basename"."sprintf("%0"length(int(scan_count/max))"d", i)"."ext
	}' $1

	echo $dir

	shift
done
