#!/usr/bin/env bash

set -e

#	(I uploaded that codes (PeptideMatch java codes) that does this. It was previously published by another group. If you like, you can modify these codes or re-write it in your language you prefer.)
#
#	                If “Modified sequence” contains S, T, or Y, then please look for the location of either S, T, or Y in protein and record it. If there are more than one of these, then record all separated by comma.
#
#	                If S, T, or Y are modified in “Modified sequence”, then please record its location in STY modification position column.
#
#	                If “Modified sequence” contains M, then please look for the location of M in protein and record it. If there are more than one site, then record them separated by comma.
#
#	                If M modified in “Modified sequence”, then please record its location in M modification position column.
#
#	 
#
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
#	 
#
#	Cartoon Sample
#
#	Fasta file:
#
#	>PROTEIN1
#
#	AAASSSTTYAAAM
#
#	 
#
#	Evidence file:
#
#	Modified Sequence Leading razor protein
#
#	_ASSST[80]TY_
#
#	_SSS[80]T[80]_
#
#	 
#
#	Then,
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
#	…
#
#	 
#
#	Thank you!!
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




#	Select only those matching the requested amino acids.

#	[[ -f evidence.STYM.txt ]] || tail -n +2 evidence.txt | awk -F"\t" '( $4 ~ /[STYM]/ ){print $4,$15}' | sort -u > evidence.STYM.txt

#	for each modified sequence;
#		remember its protein
#		find position of modified sequence and protein in protein database
#		find given amino acid positions
#		find given modified amino acid positions
#			

#	This will likely be better done in something like ruby.






script=$( basename $0 )

#	Defaults:

function usage(){
	echo
#	echo "Wrapper around calling Drop-seq_alignment.sh"
	echo
	echo "Usage:"
	echo
	echo "$script <OPTIONS>"
	echo
	echo "Options:"
	echo "  --amino_acids COMMA SEPARATED TEXT"
	echo "  --evidence TEXT TSV FILE"
	echo "  --protein FASTA FILE"
#	echo "	--estimated-num-cells (-n) INTEGER : "
#	echo "	--genomedir (-g) STRING : Directory of STAR genome directory"
#	echo "	--referencefasta (-r) STRING : Reference fasta of the Drop-seq reference metadata bundle"
#	echo "	--max INTEGER ......... Maximum number of scans per file"
#	echo "	--verbose .................. NO VALUE, just boolean flag"
#	echo
#	echo "Default option values:"
#	echo "	--estimated-num-cells ... ${num_cells}"
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
		-a|--amino_acids)
			shift; amino_acids=$1; shift;;
		-e|--evidence)
			shift; evidence=$1; shift;;
		-p|--protein)
			shift; protein_fasta=$1; shift;;
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
[[ -z $amino_acids ]] && ( echo -e "\nAmino Acids Required" && usage )
[[ -f $evidence ]] || ( echo -e "\nEvidence file Required" && usage )
[[ -f $protein_fasta ]] || ( echo -e "\nProtein fasta file Required" && usage )

protein_base_name=${protein_fasta%%.*}
echo $protein_base_name

if [[ ! -d $protein_base_name ]] ; then
	echo "Creating protein database."
	java -jar PeptideMatchCMD_1.0.jar -a index -d ${protein_fasta} -i ${protein_base_name}
fi


acids=$(echo $amino_acids | sed 's/[^[:alpha:]]//g' )
echo $acids

[[ -f evidence.${acids}.txt ]] || tail -n +2 evidence.txt | awk -F"\t" '( $4 ~ /['$acids']/ ){print $4"\t"$15}' | sort -u > evidence.${acids}.txt



echo -e -n "Modified sequence\tLeading razor protein\tStart position\tEnd position" > MatchedModification.txt
for aa in $( echo ${amino_acids} | awk -F, '{ for(i=1;i<=NF;i++) print $i }' ) ; do
	echo -e -n "\t${aa} position\t${aa} modification position" >> MatchedModification.txt
done
echo >> MatchedModification.txt


while read line; do
#	i=$[ $i + 1 ]
#	if [ $i -gt 30 ] ; then break ;fi

  echo $line
	sequence=${line%	*}
	echo $sequence
	protein=${line#*	}
	echo $protein

	cleaned_sequence=$( echo ${sequence} | sed -e 's/_//g' -e 's/([[:alpha:]]*)//g' )
	echo $cleaned_sequence

	out=${protein_base_name}.${cleaned_sequence}.txt

	[[ -f ${out} ]] ||  java -jar PeptideMatchCMD_1.0.jar -a query -i ${protein_base_name} -q ${cleaned_sequence} -o ${out}

	number_matches_in_protein=$( cat ${out} | grep ${protein} | wc -l )
	#	I suppose its possible that sequence is repeated in a given protein.
	#	Not sure what to do if it is.

	[[ ${number_matches_in_protein} -eq 0 ]] && ( echo "None of the matches match the xpected protein."; exit )

	match_start=$( cat ${out} | grep ${protein} | awk '{print $4}' )
	match_end=$( cat ${out} | grep ${protein} | awk '{print $5}' )


	echo -e -n "${sequence}\t${protein}\t${match_start}\t${match_end}" >> MatchedModification.txt

	for aa in $( echo ${amino_acids} | awk -F, '{ for(i=1;i<=NF;i++) print $i }' ) ; do
		positions=$( echo ${cleaned_sequence} | grep -o . | grep -n "[${aa}]" | awk -F: -v s=${match_start} '{x=x""$1+s";"}END{printf(substr(x, 1, length(x)-1))}' )

		all_modified_positions=$( echo ${sequence} | sed 's/_//g' | awk -v s=${match_start} '{ while(( m = match($0,/\(/) ) > 0 ){ x=x""s+m-1";"; sub(/\(.{2,4}\)/,"",$0); } }END{printf(substr(x, 1, length(x)-1)) }' )

		
		#	comm -12 <( echo "12;14;21" | sed 's/;/\n/g' | sort ) <( echo "2;14"| sed 's/;/\n/g' | sort )
		#	comm is meant to compare 2 files containing sorted lists.
		#		-12 means ignore items only in 1 of those lists, leaving just those that are in both.

		modified_positions=$( comm -12 <( echo ${positions} | sed 's/;/\n/g' | sort ) <( echo ${all_modified_positions} | sed 's/;/\n/g' | sort ) )


		echo -e -n "\t${positions}\t${modified_positions}" >> MatchedModification.txt

	done

	echo >> MatchedModification.txt




done < evidence.${acids}.txt
#done < <( head evidence.${acids}.txt )

echo 'Done'

