#!/usr/bin/env bash

set -e

#	The second output file (name it as ProteinModification.txt) contains 1) Leading razor protein, 2) AA, 3) AA location, 4) Sequences with Modified AA, and 5) Sequences with Unmodified AA
#	For each leading razor protein that were observed in at least once in evidence file, find all S, T, Y, or M. AA will contain one of these amino acids located in a specific location.
#
#	Each location of S, T, Y, or M will be one line. “AA location” will contain the location of the corresponding amino acid (one location per line).
#
#	For “Sequences with Modified AA”, records all unique sequences with the corresponding modified “AA”.
#
#	For “Sequences with Unmodified AA”, records all unique sequences with the corresponding unmodified “AA”.
#
#	 
#
#	Cartoon Sample
#
#		Fasta file:
#
#		>PROTEIN1
#		AAASSSTTYAAAM
#
#
#	Evidence file:
#
#		Modified Sequence 	Leading razor protein
#		_ASSST[80]TY_			yadayadayada
#		_SSS[80]T[80]_		somethingelse
#
#
#	Output file 1:
#
#	Modified sequence\t Leading razor protein\t Start position\t End position\t STY position\t STY modification position\t M position\t M modification position \n
#
#	_ASSST[80]TY_\t       PROTEIN1\t     3\t      9\t    4;5;6;7;8;9\t    7\t        \t      \t
#
#	_SSS[80]T[80]_\t      PROTEIN1\t     4\t      7\t    4;5;6;7\t        6;7\t      \t      \t
#
#	 
#
#	Output file 2:
#
#	Leading razor protein\t AA\t AA location\t Sequences with Modified AA\t Sequences with Unmodified AA\t
#
#	PROTEIN1\t               S\t           4\t     \t                       _ASSST[80]TY_;_SSS[80]T[80]_
#	PROTEIN1\t               S\t           5\t     \t                       _ASSST[80]TY_;_SSS[80]T[80]_
#	PROTEIN1\t               S\t           6\t     _SSS[80]T[80]_\t         _ASSST[80]TY_
#	PROTEIN1\t               T\t           7\t\    _ASSST[80]TY_; _SSS[80]T[80]_\t
#
#
#	I'll see if I can just uniq the modified sequence / leading razor protein from within the script itself.
#
#	All positions are relative to the start of the full protein found in the fasta file.
#
#	Ignore the underscores.
#
#	Modifications come after the column.
#
#	The modification value itself is irrelevant in this script.
#
#
#	modification_locator.bash --amino_acids STY,M --evidence evidence.txt --protein uniprot-organism+homo+sapiens.fasta
#
#	As each record in ProteinModification.txt, comes from multiple records in the source,
#	it will likely be simplest to just make it from MatchModification.txt.
#
#	protein_modification_locator.bash --match MatchedModification.txt








script=$( basename $0 )

#	Defaults:
match_file="MatchedModification.txt"

function usage(){
	echo
#	echo "Wrapper around calling Drop-seq_alignment.sh"
	echo
	echo "Usage:"
	echo
	echo "$script <OPTIONS>"
	echo
	echo "Options:"
	echo "	--match TSV TEXT FILE"
#	echo "	--estimated-num-cells (-n) INTEGER : "
#	echo "	--genomedir (-g) STRING : Directory of STAR genome directory"
#	echo "	--referencefasta (-r) STRING : Reference fasta of the Drop-seq reference metadata bundle"
#	echo "	--max INTEGER ......... Maximum number of scans per file"
#	echo "	--verbose .................. NO VALUE, just boolean flag"
#	echo
	echo "Default option values:"
	echo "	--match ... ${match_file}"
#	echo "	--genomedir ............. ${genomedir}"
#	echo "	--referencefasta ........ ${referencefasta}"
#	echo "	--max ........ ${max}"
#	echo "	--verbose .......... ${verbose}"
	echo
	echo "Examples:"
#	echo "	$script file.bam"
#	echo "	$script --max 10 file.mgf"
	echo
	echo
	exit 1
}

while [ $# -ne 0 ] ; do
	case $1 in
		-m|--match)
			shift; match_file=$1; shift;;
#		-e|--evidence)
#			shift; evidence=$1; shift;;
#		-p|--protein)
#			shift; protein_fasta=$1; shift;;
#		-n|--estimated-num-cells)
#			shift; num_cells=$1; shift;;
#		-g|--genomedir)
#			shift; genomedir=$1; shift;;
#		-r|--referencefasta)
#			shift; referencefasta=$1; shift;;
#		-m|--m*)
#			shift; max=$1; shift;;
#		-v|--v*)
#			verbose=true; shift;;
		-*)
			echo ; echo "Unexpected args from: ${*}"; usage ;;
		*)
			break;;
	esac
done

#	Basically, this is TRUE AND DO ...
#[ $# -eq 0 ] && usage
if [[ ! -f $match_file ]] ; then
	echo -e "\nMatchModification file Required"
	usage
fi

date

echo -e "Leading razor protein\tAA\tAA location\tSequences with Modified AA\tSequences with Unmodified AA" > ProteinModification.txt

#	Find out which amino acids and groups are in match file
#amino_acids="STY,M"
amino_acids=$(head -1 ${match_file} | tr "\t" "\n" | grep "modification position" | awk '{print $1}' | tr "\n" ",")
amino_acids=${amino_acids%,}
echo $amino_acids

proteins=$( tail -n +2 ${match_file} | awk -F"\t" '{print $2}' | sort -u )

for protein in ${proteins} ; do

	matches=$( awk -F"\t" '( $2 == "'${protein}'"){ print $0 }' ${match_file} )

#echo $matches | awk -F"\t" '{print 





	echo ${protein} >> ProteinModification.txt




done

echo 'Done'
date
