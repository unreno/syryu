#!/usr/bin/env bash

script=`basename $0`
script_dir=`dirname $0`

#	Defaults:

function usage(){
	echo
	echo "Outputs file of same name with .SCANS.mgf appended to root name."
	echo
	echo "reads in an MGF file"
	echo "outputs input MGF file with added SCANS=# line parsed from TITLE."
	echo
	echo "Usage:"
	echo
	echo "$script <OPTIONS> mgf_file(s)"
	echo
	echo "Options:"
#	echo "	--verbose .................. NO VALUE, just boolean flag."
#	echo "	          .................. Output file will be invalid MGF as contains a lot of other stuff"
	echo
	echo "Default option values:"
#	echo "	--verbose ........... ${verbose}"
	echo
	echo "Notes:"
	echo
	echo "Examples:"
	echo "	$script 140521_EOC_MCis_NT_3.mgf"
	echo
	echo
	exit 1
}

while [ $# -ne 0 ] ; do
	case $1 in
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


#	-v verbose=$verbose \
awk \
'
BEGIN{ }
( /^BEGIN IONS/ ){ scans="" }
( /^TITLE=/ ){ 
	match($0,/scan=[[:digit:]]+/)
	if( RSTART !=0 ) scans=substr($0,RSTART+5,RLENGTH-5)
}
( /^CHARGE=/ ){ print "SCANS="scans }
{	print }' $1 > ${1%.*}.SCANS.mgf

	shift

done


