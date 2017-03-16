#!/usr/bin/env bash


while [ $# -ne 0 ] ; do
#	echo $1

awk '
function split_and_print() {
	print name
	print sequence
}
( !/^>/ ){
	sequence=sequence$0
}
( /^>/ ){
	if( name ) split_and_print()
	
	name=$1
	gsub(/^>/,"",name)
	sequence=""
}
END{
	split_and_print()
}' $1

	shift
done

