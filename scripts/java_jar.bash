#!/usr/bin/env bash

script=`basename $0`

#	Defaults:
verbose=false

function usage(){
	echo
	echo "Calls java and searches CLASSPATH for jar file and adds path prefix if found."
	echo
	echo "Usage:"
	echo
	echo "$script <OPTIONS> just_one_jar_file"
	echo
	echo "Options:"
	echo "	--verbose .................. NO VALUE, just boolean flag"
	echo
	echo "Default option values:"
	echo "	--verbose .......... ${verbose}"
	echo
	echo "Examples:"
	echo "	$script "
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
#[ $# -eq 0 ] && usage
[ $# -ne 1 ] && usage





#
jar_file=$1


#	set CLASSPATH=D:\myprogram;D:\myprogram\lib\supportLib.jar
#	split on ; and rejoin with " "
#	add "." so searches current dir too?
#	perhaps $PATH?

#	search for jar file in CLASSPATH
#find ~/igv ~/tmpmount ~/gdc_api/ -name *txt

[ -z $CLASSPATH ] && echo "Your CLASSPATH is blank. Searching from current working dir."

: ${CLASSPATH:="."}


#	${CLASSPATH//;/ } replaces ALL ;'s with ' 's
#	2>/dev/null so no "Permission denied" lines

jar_file_with_path=$( find ${CLASSPATH//;/ } -name "${jar_file}" 2>/dev/null | head -1 )


cmd="java -jar \"$jar_file_with_path\""

if $verbose ; then
	echo "Calling ${cmd}"
fi

eval $cmd

