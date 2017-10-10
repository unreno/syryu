#!/usr/bin/env ruby

require 'csv'
require 'optparse'

def usage
	puts
	puts "Usage"
	puts
	puts
	puts
	puts
	puts "#{$0} <formatted moda output file> <observed spectra mgf file>"
	puts
	exit
end


options = { suffix: 'TESTING'}
OptionParser.new do |opt|
#	opt.on('--first_name FIRSTNAME') { |o| options[:first_name] = o }
#	opt.on('--last_name LASTNAME') { |o| options[:last_name] = o }
	opt.on('--suffix TESTING3') { |o| options[:suffix] = o }
end.parse!

#puts options

usage if ARGV.length < 2




#date=$(date "+%Y%m%d%H%M%S")
#date="TESTING"
date=options[:suffix]
theodir="#{ARGV[0]}_theospec_#{date}/"
Dir.mkdir(theodir) unless Dir.exists?(theodir)
#obsdir="#{ARGV[0]}_obsspec_#{date}/"
#Dir.mkdir(obsdir) unless Dir.exists?(obsdir)


marshaled_observed_spectra="#{ARGV[1]}.#{date}.marshal"
if File.exists?( marshaled_observed_spectra )
	observed_spectra=Marshal.load(File.binread(marshaled_observed_spectra))
else
	observed_spectra=[]
	File.open( ARGV[1], 'r' ).each do |line|
	#	@spectrum={} if @spectrum.nil?
		if line =~ /BEGIN IONS/
			@spectrum={}
			#	nothing
		elsif line =~ /^TITLE=(.*)$/
			@spectrum['title'] = $1
		elsif line =~ /^RTINSECONDS=(.*)$/
			@spectrum['rtinseconds'] = $1.to_f
		elsif line =~ /^PEPMASS=(.*)$/
			@spectrum['pepmass'] = $1
		elsif line =~ /^SCANS=(.*)$/
			@spectrum['scans'] = $1.to_i
		elsif line =~ /^CHARGE=(.*)$/
			@spectrum['charge'] = $1
		elsif line =~ /END IONS/
			observed_spectra.push @spectrum
			@spectrum = nil
		elsif line =~ /^(.+)\s+(.+)$/
			@spectrum['mzis'] = [] unless @spectrum['mzis']
			@spectrum['mzis'].push([$1.to_f,$2.to_f])
		end
	end
	File.open(marshaled_observed_spectra, 'wb') {|f| f.write(Marshal.dump(observed_spectra))}
end



marshaled_formatted_moda_output="#{ARGV[0]}.#{date}.marshal"
if File.exists?( marshaled_formatted_moda_output )
	formatted_moda_output=Marshal.load(File.binread(marshaled_formatted_moda_output))
else
	formatted_moda_output = CSV.read( ARGV[0],'rb',{ col_sep: "\t", headers: true })
	formatted_moda_output.each do |row|
	
	#<CSV::Row "output_index":"2" "spec_index":"19" "observed_MW":"3094.4649" "charge_state":"5" "scan_number":"1535" "rank":"1" "calculated_MW":"3094.4650" "delta_mass":"-0.0001" "score":"46" "probability":"0.3989" "peptide":"R.KTNDKDEKKEDGKQAENDSSNDDKTKK.S" "protein":"sp|Q9BXP5" "pept_position":"301~327" "mod1":"NA" "mod2":"NA" "mod3":"NA" "PlainPeptide":"KTNDKDEKKEDGKQAENDSSNDDKTKK">
	
		#	as peptide could include decimals,
		#	split on decimals to pick off the first and last and then parse the middle
		parts = row['peptide'].split(/\./); #	-.QWLPKD-123EESFLQFKYQALQVP+117.P	(theoretically could include decimals)
		a=parts[0]               #	-
		b=parts[1..-2].join('.') #	QWLPKD-123EESFLQFKYQALQVP+117
		c=parts[-1]              #	P
	
		#	parse middle extracting numbers and replacing them with codes passed to theospec
		numbers=b.split(/[A-Z]+/)    #	["",-123,+117]
		letters=b.split(/[0-9.+-]+/) #	[ "QWLPKD","EESFLQFKYQALQVP"]
	
		outpeptide=""
	
		letters.each_index do |i|
			outpeptide="#{outpeptide}#{letters[i]}"
			outpeptide="#{outpeptide}#{i}" unless( numbers[i+1].to_s.empty? )
		end
	
		#	"QWLPKD1EESFLQFKYQALQVP2"
	
		z=( row['charge_state'].to_i > 5 ) ? 5 : row['charge_state']
	
		theospec_command="theospec -z#{z}"
		numbers.delete_if{|n|n.empty?}.each_with_index do |number,i|
			theospec_command="#{theospec_command} -c#{i}=#{number}"
		end
		#	-c0=-123 -c1=+117
	
		theospec_output="#{theodir}#{row['scan_number']}.#{row['rank']}.theospec.txt"
	
		theospec_command="#{theospec_command} #{a}.#{outpeptide}.#{c} > #{theospec_output}"
	#	puts theospec_command
		unless File.exists?( theospec_output )
	#		puts "Running theospec"
			system(theospec_command)
		end
	
#		row['theospec']=File.readlines( theospec_output ).collect{|l|l.chomp}
		row['theospec_comments']=""
		row['theospec']=[]
		File.open( theospec_output, 'r' ).each do |line|
			next if line.empty?
			if line =~ /^[[:digit:]]/
				#	Not sure what these all are
				#158.092403080(10);+1;y1;0;(-NH3)	
				parts=line.split(';')
				row['theospec'].push({
					mz: parts[0].to_f,  #	158.092403080(10) to_f -> 158.092403080
					b: parts[1],        #	+1
					ion: parts[2],      #	y1
					d: parts[3],        #	0
					e: parts[4].rstrip, #	(-NH3)
				})
			else
				row['theospec_comments'] << "#{line}\n"
			end
		end
	
	#	row['theospec']=`#{theospec_command} #{a}.#{outpeptide}.#{c}`
	#	puts row['theospec'].inspect
	#	puts theospec_command
	
	end
	
	File.open(marshaled_formatted_moda_output, 'wb') {|f| f.write(Marshal.dump(formatted_moda_output))}
end



csvout = CSV.open( "#{ARGV[0]}.#{date}.ion_statistics",'w',{ col_sep: "\t" })

csvout.puts %w{scan_number rank number_observed_mz number_theoretical_mz max_intensity nb ny b_i intensity_b_i y_i intensity_y_i b_error y_error}

formatted_moda_output.each do |row|
	#<CSV::Row "output_index":"2" "spec_index":"19" "observed_MW":"3094.4649" "charge_state":"5" "scan_number":"1535" "rank":"1" "calculated_MW":"3094.4650" "delta_mass":"-0.0001" "score":"46" "probability":"0.3989" "peptide":"R.KTNDKDEKKEDGKQAENDSSNDDKTKK.S" "protein":"sp|Q9BXP5" "pept_position":"301~327" "mod1":"NA" "mod2":"NA" "mod3":"NA" "PlainPeptide":"KTNDKDEKKEDGKQAENDSSNDDKTKK">

	puts "#{row['scan_number']} : #{row['rank']}"
	observed = observed_spectra.find{|s| s['scans'] == row['scan_number'].to_i }		#	FIRST 
	number_observed_mz = observed['mzis'].length
	puts "The total number of observed m/z values : #{number_observed_mz}"
	number_theoretical_mz = row['theospec'].length
	puts "The total number of theoretical m/z values : #{number_theoretical_mz}"
	max_intensity=observed['mzis'].collect{|o|o[1]}.max
	puts "The largest intensity in spectrum whether it is matched or not : #{max_intensity}"

	puts "Matching ..."
	omt = observed['mzis'].collect do |o|
#		puts o[0]
		#	matching algorithm likely needs fine tuning
#	No duplicate matches. Remove match as matched
#	use find_index and the delete_at , which returns the now deleted value

#	Try find_all _indexes? sort by difference and intensity??????
#	.each_index.select{|i| ...[i] == ??? } => an array of indexes
		
		index = row['theospec'].find_index{|r| ( r[:mz] > ( o[0] - 0.5 ) ) && ( r[:mz] < ( o[0] + 0.5 ) ) }
		matched = ( index ) ? row['theospec'].delete_at(index) : {}

		{ mz: o[0], int: o[1], matched: matched }
#			matched: row['theospec'].select{|r| ( r[:a] > ( o[0] - 0.5 ) ) && ( r[:a] < ( o[0] + 0.5 ) ) } || []
#		puts matches.inspect
	end
	puts omt.inspect

	bions = omt.select{|o| o[:matched][:ion] =~ /^b/ }
	yions = omt.select{|o| o[:matched][:ion] =~ /^y/ }

	nb = bions.length
	puts "# of observed b-ion m/z values matched (b-ion is starting with b in theoretical spec; for example, b2, b3, b4, ….) Let’s call it as “Nb” : #{nb}"

	ny = yions.length
	puts "# of observed y-ion m/z values matched (y-ion is starting with y in theoretical spec; for example, y2, y5,…) Let’s call it as “Ny” : #{ny}"

	b_i = bions.collect{|i|i[:int]}.join(';')
	puts "The list of matched observed b-ion intensities in the spectrum (in mgf files). – separated by semicolon : #{b_i}"

	intensity_b_i = bions.collect{|i|i[:int]*100/max_intensity}.join(';')
	puts "The list of matched observed b-ion intensities divided by the largest intensity in spectrum, and times 100. Let’s call it “Intensity_b_i” for each ith observed peak.  – separated by semicolon : #{intensity_b_i}"

	y_i = yions.collect{|i|i[:int]}.join(';')
	puts "The list of matched observed y-ion intensities (in mgf files) – separated by semicolon : #{y_i}"

	intensity_y_i = yions.collect{|i|i[:int]*100/max_intensity}.join(';')
	puts "The list of matched observed y-ion intensities divided by the largest intensity in spectrum, and times 100. Let’s call it “Intensity_y_i” for each ith observed peak. – separated by semicolon : #{yions.collect{|i|i[:int]*100/max_intensity}.join(';')}"



#	m/z error = observed m/z minus theoretical m/z

	b_error = 0
#-m/z errors between observed and theoretical m/z values for b-ion matches – separated by semicolon

	y_error = 0
#-m/z errors between observed and theoretical m/z values for y-ion matches – separated by semicolon



	csvout.puts [row['scan_number'], row['rank'], number_observed_mz, number_theoretical_mz, max_intensity, nb, ny, b_i, intensity_b_i, y_i, intensity_y_i,b_error,y_error]

	puts
	puts
	puts
	puts
end

csvout.close

