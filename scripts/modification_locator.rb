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
	puts "Produces MatchedModification.txt, ProteinModification.txt, "
	puts "MatchedModificationNONE.txt and MatchedModificationMULTIPLE.txt"
	puts
	puts "#{$0} [options] "#<formatted moda output file> <observed spectra mgf file>"
	puts
	puts "Options:"
	puts "	--amino_acids COMMA SEPARATED TEXT"
	puts "	--evidence TEXT TSV FILE"
	puts "	--protein FASTA FILE"
#	puts "	--suffix STRING ...... working directory suffixes"
#	puts "	--log_count INTEGER .. number of records logged before silencing output"
	puts
	puts "Defaults/current values:"
#	puts "	--suffix #{options[:suffix]}"
	puts
	puts "Examples:"
#	puts "#{$0} out_blind_140521_EOC_MCis_T2_3.txt.OUTPUT.tsv 140521_EOC_MCis_T2_3.SCANS.mgf"
	puts "#{$0} --amino_acids STY,M --evidence evidence.txt --protein uniprot-organism+homo+sapiens.fasta"
	puts
	exit
end


######################################################################


#	Setup future redirect of STDOUT if desired
#stdout_redirected = false
#original_stdout = $stdout.clone
#	In the future, could redirect STDOUT to /dev/null
#	if !stdout_redirected && options.has_key?(:log_count) && record_number > options[:log_count].to_i
#		puts "Record number(#{record_number}) has exceeded requested log count(#{options[:log_count]})."
#		puts "Redirecting the rest of the output to /dev/null."
#		stdout_redirected = true
#		$stdout.reopen("/dev/null", "a")
#	end


#	Must be called before option parsing as they remove the items.
puts "Command: #{$0} #{$*.join(' ')}"


options = {}
#options = { suffix: 'TESTING'}
#options = { suffix: Time.now.strftime("%Y%m%d%H%M%S") }
OptionParser.new do |opt|
	opt.on("--amino_acids COMMA SEPARATED TEXT"){|o| options[:amino_acids] = o }
	opt.on("--evidence TEXT TSV FILE"){|o| options[:evidence_file] = o }
	opt.on("--protein FASTA FILE"){|o| options[:protein_file] = o }
#	opt.on('--suffix TESTING') { |o| options[:suffix] = o }
#	opt.on('--log_count 10') { |o| options[:log_count] = o }
end.parse!

puts options


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

class String

	def indices_of_chars(chars)
    start, indices = -1, []
    indices << start while start = (self.index /[#{chars}]/, start + 1)
    indices
	end

	def indices_of_mods()
		tmp=self
    indices = []
		while (i = tmp.index /\(/)
    	indices << i - 1
			tmp.sub! /\(.*?\)/,''
		end
    indices
	end

end



class Evidence

	attr_accessor :sequence, :protein, :cleaned_sequence, :alignments, :modified_positions

	@@dir = "PeptideMatchOutput"

	def initialize(h={})
		@sequence=h["Modified sequence"]
		@protein=h["Leading razor protein"]
		@cleaned_sequence=sequence.gsub(/_/,'').gsub(/\([[:alpha:]]*\)/,'')
		@alignments = []
		@modified_positions=sequence.gsub(/_/,'').indices_of_mods
	end

	def query( database )
		Dir.mkdir(@@dir) unless Dir.exists?(@@dir)
		outfile="#{@@dir}/#{database}.#{cleaned_sequence}.txt"

		#	puts outfile

		if !File.exists?( outfile )
			puts "Querying."
			system("java -jar PeptideMatchCMD_1.0.jar -a query -i #{database} -q #{cleaned_sequence} -o #{outfile}")
		else
			puts "Catting existing query."
		end
		
		File.open(outfile,'rb').each do |line|
			next if line =~ /^#/
			l = line.split
			if l[0] == cleaned_sequence && l[1] == protein
				##Query	Subject	SubjectLength	MatchStart	MatchEnd
				alignments.push({ match_start: l[3].to_i, match_end: l[4].to_i,
					absolute_positions: {}, modified_positions: [] })
			end
		end
	end

end

##################################################

protein_base_name=options[:protein_file].gsub(/.[^.]*$/,'').gsub(/.*\/(.*)$/,'\1')

if !Dir.exists?( protein_base_name )
	puts "Creating #{protein_base_name} protein database."
	system("java -jar PeptideMatchCMD_1.0.jar -a index -d #{options[:protein_file]} -i #{protein_base_name}")
else
	puts "#{protein_base_name} protein database exists."
end

##################################################

acids = options[:amino_acids].gsub(/[^[:alpha:]]/,'')
select_evidence_file = options[:evidence_file].gsub(/(.[^.]+)$/,".#{acids}\\1")
if !File.exists?( select_evidence_file )
	puts "Filtering evidence."
	acids_regex=Regexp.new("[#{acids}]")	
	select_evidence = CSV.open(options[:evidence_file],'rb',
			{ headers: true, col_sep: "\t"  }).collect do |line|
		[line["Modified sequence"],line["Leading razor protein"]] if line["Modified sequence"].match(acids_regex)
	end.compact

	select_evidence_csv = CSV.open(select_evidence_file,'w', { col_sep: "\t"  })
	select_evidence_csv.puts ["Modified sequence","Leading razor protein"]
	select_evidence.sort.uniq.each do |line|
		select_evidence_csv.puts line
	end
	select_evidence_csv.close
else
	puts "Evidence already filtered."
end

##################################################

puts "Starting"
puts Time.now

evidences=[]

no_matches = CSV.open("MatchedModificationNONE.txt",'w', {col_sep: "\t" })
no_matches.puts ["Modified sequence","Leading razor protein"]
multiple_matches = CSV.open("MatchedModificationMULTIPLE.txt",'w', {col_sep: "\t" })
multiple_matches.puts ["Modified sequence","Leading razor protein"]
matched_mod = CSV.open("MatchedModification.txt",'w', {col_sep: "\t" })
matched_mod_header = ["Modified sequence","Leading razor protein","Start position","End position"]
options[:amino_acids].split(",").each do |acid|
	matched_mod_header.push "#{acid} position"
	matched_mod_header.push "#{acid} modification position"
end
matched_mod.puts matched_mod_header

c = CSV.open(select_evidence_file,'rb')
record_count = c.readlines.size
c.close
(c=CSV.open(select_evidence_file,'rb',
	{ headers: true, col_sep: "\t"  })).each do |line|

	puts "Processing #{c.lineno}/#{record_count} - #{line['Modified sequence']}:#{line['Leading razor protein']}"

	e = Evidence.new(line.to_hash)
	e.query( protein_base_name )
	evidences.push e

	if e.alignments.length < 1
		puts "None of the matches match the expected protein."
		no_matches.puts line
	else
		if e.alignments.length > 1
			puts "More than one matches match the expected protein. Reporting all."
			multiple_matches.puts line
		end

		e.alignments.each do |alignment|

			matched_mod_line = [e.sequence,e.protein,alignment[:match_start],alignment[:match_end]]
	
			acids.split(//).each do |acid|
				positions = e.cleaned_sequence.indices_of_chars(acid).collect{|i| alignment[:match_start] + i }
				positions.each { |position| alignment[:absolute_positions][position] = acid }
			end

			options[:amino_acids].split(",").each do |acid_group|

				absolute_positions=alignment[:absolute_positions].select{|k,v| acid_group.match(v) }.collect{|k,v| k }

				modified_positions=e.modified_positions.collect{|i| alignment[:match_start] + i } & absolute_positions
				alignment[:modified_positions] += modified_positions

				matched_mod_line << absolute_positions.join(";")
				matched_mod_line << modified_positions.join(";")
			end

			matched_mod.puts matched_mod_line

		end	#	e.alignments.each do |alignment|

	end	#	e.alignments.length >= 1

	#	Can't use $. as it is the last counter used and in this case, its the query output, not the csv file.
	#break if $. > 10
#	break if c.lineno > 10

end		#	CSV.open(options[:evidence_file],'rb',

no_matches.close
multiple_matches.close
matched_mod.close

##################################################

protein_mod=CSV.open("ProteinModification.txt",'w', {col_sep: "\t" })
protein_mod.puts ["Leading razor protein","AA","AA location","Sequences with Modified AA","Sequences with Unmodified AA"]


#	Leading razor protein\t AA\t AA location\t Sequences with Modified AA\t Sequences with Unmodified AA\t
#
#	PROTEIN1\t               S\t           4\t     \t                       _ASSST[80]TY_;_SSS[80]T[80]_
#	PROTEIN1\t               S\t           5\t     \t                       _ASSST[80]TY_;_SSS[80]T[80]_
#	PROTEIN1\t               S\t           6\t     _SSS[80]T[80]_\t         _ASSST[80]TY_
#	PROTEIN1\t               T\t           7\t\    _ASSST[80]TY_; _SSS[80]T[80]_\t


evidences.collect{|e|e.protein}.uniq.sort.each do |protein|

	evidences_for_this_protein = evidences.select{|e|e.protein == protein}

	positions = {}

	alignments = evidences_for_this_protein.collect do |e|
		e.alignments.each do |alignment| 
			positions.update alignment[:absolute_positions]
			alignment[:sequence] = e.sequence
		end
		e.alignments
	end.flatten

	positions.keys.sort.each do |position|

		acid=positions[position]

		alignments_at_this_position = alignments.select{|a| a[:match_start] <= position && a[:match_end] >= position }

		protein_mod_line = [ protein, acid, position, 
			alignments_at_this_position.select{|a|  a[:modified_positions].include? position
				}.collect{|a| a[:sequence] }.join(';'),
			alignments_at_this_position.select{|a| !a[:modified_positions].include? position
				}.collect{|a| a[:sequence] }.join(';') ]

		protein_mod.puts protein_mod_line

	end
	
end	#	evidences.collect{|e|e.protein}.uniq.sort.each do |protein|

protein_mod.close

puts "Done"
puts Time.now
