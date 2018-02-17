#!/usr/bin/env ruby

require 'csv'
require 'optparse'

def usage(options={})
	puts
	puts "Usage"
	puts
#	puts "Uses a date suffix for all created files."
#	puts "If going to run multiple times, for whatever reason, you can specify the same suffix with suffix option."
#	puts "This way, most of the parsing and the theospec calls won't need rerun."
	puts
	puts "#{$0} [options] "#<formatted moda output file> <observed spectra mgf file>"
	puts
	puts "Options:"
	puts "	--amino_acids COMMA SEPARATED TEXT"
	puts "	--evidence TEXT TSV FILE"
	puts "	--protein FASTA FILE"
#	puts "	--suffix STRING ...... output file and directory suffixes"
#	puts "	--log_count INTEGER .. number of records logged before silencing output"
	puts
	puts "Defaults/current values:"
#	puts "	--suffix #{options[:suffix]}"
	puts
#	puts "#{$0} out_blind_140521_EOC_MCis_T2_3.txt.OUTPUT.tsv 140521_EOC_MCis_T2_3.SCANS.mgf"
	puts
	exit
end


######################################################################


stdout_redirected = false
original_stdout = $stdout.clone
#	$stdout.reopen("/dev/null", "a")


#	Must be called before option parsing as they remove the items.
puts "Command: #{$0} #{$*.join(' ')}"


#options = { suffix: 'TESTING'}
options = { suffix: Time.now.strftime("%Y%m%d%H%M%S") }
OptionParser.new do |opt|
	opt.on("--amino_acids COMMA SEPARATED TEXT"){|o| options[:amino_acids] = o }
	opt.on("--evidence TEXT TSV FILE"){|o| options[:evidence_file] = o }
	opt.on("--protein FASTA FILE"){|o| options[:protein_file] = o }
#	opt.on('--suffix TESTING') { |o| options[:suffix] = o }
#	opt.on('--log_count 10') { |o| options[:log_count] = o }
end.parse!

#puts options


######################################################################

#puts "Using suffix '#{options[:suffix]}'"



if options[:amino_acids].to_s.empty?
	puts "\nAmino Acids Required"
	usage(options)
end
if options[:evidence_file].to_s.empty? or !File.exists?(options[:evidence_file])
	puts "\nEvidence file Required"
	usage(options)
end
if options[:protein_file].to_s.empty? or !File.exists?(options[:protein_file])
	puts "\nProtein fasta file Required"
	usage(options)
end



#usage(options) #if ARGV.length < 2 



exit



#matched_mod=File.open("MatchedModification.txt",'w')
#matched_mod.close


protein_mod=File.open("ProteinModification.txt",'w')
protein_mod.puts "Leading razor protein\tAA\tAA location\tSequences with Modified AA\tSequences with Unmodified AA" 

CSV.open("MatchedModification.txt",'rb',
	{ headers: true, col_sep: "\t"  }).each do |line|

	puts line.to_hash





end






protein_mod.close





#	protein_base_name=${protein_fasta%%.*}
#	echo $protein_base_name
#	
#	if [[ ! -d $protein_base_name ]] ; then
#		echo "Creating protein database."
#		java -jar PeptideMatchCMD_1.0.jar -a index -d ${protein_fasta} -i ${protein_base_name}
#	fi
#	
#	
#	acids=$(echo $amino_acids | sed 's/[^[:alpha:]]//g' )
#	echo $acids
#	
#	[[ -f evidence.${acids}.txt ]] || tail -n +2 evidence.txt | awk -F"\t" '( $4 ~ /['$acids']/ ){print $4"\t"$15}' | sort -u > evidence.${acids}.txt
#	
#	
#	
#	echo -e -n "Modified sequence\tLeading razor protein\tStart position\tEnd position" > MatchedModification.txt
#	for aa in $( echo ${amino_acids} | awk -F, '{ for(i=1;i<=NF;i++) print $i }' ) ; do
#		echo -e -n "\t${aa} position\t${aa} modification position" >> MatchedModification.txt
#	done
#	echo >> MatchedModification.txt
#	
#	
#	while read line; do
#	
#	  echo $line
#		sequence=${line%	*}
#		echo $sequence
#		protein=${line#*	}
#		echo $protein
#	
#		cleaned_sequence=$( echo ${sequence} | sed -e 's/_//g' -e 's/([[:alpha:]]*)//g' )
#		echo $cleaned_sequence
#	
#		out=${protein_base_name}.${cleaned_sequence}.txt
#	
#		[[ -f ${out} ]] ||  java -jar PeptideMatchCMD_1.0.jar -a query -i ${protein_base_name} -q ${cleaned_sequence} -o ${out}
#	
#		number_matches_in_protein=$( cat ${out} | grep ${protein} | wc -l )
#	
#		echo "Found ${number_matches_in_protein}"
#	
#		if [[ ${number_matches_in_protein} -eq 0 ]] ; then
#			echo "None of the matches match the expected protein."
#			echo $line >> MatchedModificationNONE.txt
#			continue
#		fi
#	
#		if [[ ${number_matches_in_protein} -gt 1 ]] ; then
#			echo "More than one matches match the expected protein."
#			echo $line >> MatchedModificationMULTIPLE.txt
#			continue
#		fi
#	
#		match_start=$( cat ${out} | grep ${protein} | awk '{print $4}' )
#		match_end=$( cat ${out} | grep ${protein} | awk '{print $5}' )
#	
#		echo -e -n "${sequence}\t${protein}\t${match_start}\t${match_end}" >> MatchedModification.txt
#	
#		for aa in $( echo ${amino_acids} | awk -F, '{ for(i=1;i<=NF;i++) print $i }' ) ; do
#			positions=$( echo ${cleaned_sequence} | fold -w1 | grep -n "[${aa}]" | awk -F: -v s=${match_start} '{printf($1+s";")}' )
#			positions=${positions%;}
#	
#			all_modified_positions=$( echo ${sequence} | sed 's/_//g' | awk -v s=${match_start} '{ while(( m = match($0,/\(/) ) > 0 ){ printf(s+m-1";"); sub(/\(.{2,4}\)/,"",$0); } }' )
#			all_modified_positions=${all_modified_positions%;}
#	
#			
#			modified_positions=$( comm -12 <( echo ${positions} | sed 's/;/\n/g' | sort ) <( echo ${all_modified_positions} | sed 's/;/\n/g' | sort ) | tr "\n" ";" )
#	
#			modified_positions=${modified_positions%;}	#	remove last semicolon
#	
#			echo -e -n "\t${positions}\t${modified_positions}" >> MatchedModification.txt
#	
#		done
#	
#		echo >> MatchedModification.txt
#	
#	
#	
#	
#	done < evidence.${acids}.txt
#	
#	echo 'Done'
#	date
#	
