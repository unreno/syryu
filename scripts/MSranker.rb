#!/usr/bin/env ruby

require 'csv'
require 'optparse'

def usage
	puts
	puts "Usage"
	puts
	puts "Uses a date suffix for all created files."
	puts "If going to run multiple times, for whatever reason, you can specify the same suffix with suffix option."
	puts "This way, most of the parsing and the theospec calls won't need rerun."
	puts
	puts "#{$0} <formatted moda output file> <observed spectra mgf file>"
	puts
	puts "MSranker.rb out_blind_140521_EOC_MCis_T2_3.txt.OUTPUT.tsv 140521_EOC_MCis_T2_3.SCANS.mgf"
	puts
	exit
end


#options = { suffix: 'TESTING'}
options = { suffix: Time.now.strftime("%Y%m%d%H%M%S") }
OptionParser.new do |opt|
#	opt.on('--first_name FIRSTNAME') { |o| options[:first_name] = o }
#	opt.on('--last_name LASTNAME') { |o| options[:last_name] = o }
	opt.on('--suffix TESTING') { |o| options[:suffix] = o }
end.parse!

#puts options
puts "Using suffix '#{options[:suffix]}'"

usage if ARGV.length < 2

#moda_base = File.basename( ARGV[0], ".*" )	#	loses path
moda_base = ARGV[0].sub(/#{File.extname(ARGV[0])}$/,'')
#mgf_base = File.basename( ARGV[1], ".*" )	#	loses path
mgf_base = ARGV[1].sub(/#{File.extname(ARGV[1])}$/,'')

theodir="#{ARGV[0]}_theospec_#{options[:suffix]}/"
Dir.mkdir(theodir) unless Dir.exists?(theodir)


#	Let’s set the default m/z error as tol=0.5. (this is used when matching observed to theoretical)
default_tol = 0.5


##################################################

#	Read, parse and store observed spectra.

marshaled_observed_spectra="#{ARGV[1]}.#{options[:suffix]}.marshal"
if File.exists?( marshaled_observed_spectra )
	observed_spectra=Marshal.load(File.binread(marshaled_observed_spectra))
else
	observed_spectra=[]
	File.open( ARGV[1], 'r' ).each do |line|
		if line =~ /BEGIN IONS/
			@spectrum={}
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

##################################################

#	Read and parse moda output, then compute and store theoretical spectra.

marshaled_formatted_moda_output="#{ARGV[0]}.#{options[:suffix]}.marshal"
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
		system(theospec_command) unless File.exists?( theospec_output )

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
					charge_state: parts[1],        #	+1
					ion: parts[2],      #	y1
					d: parts[3],        #	0
					e: parts[4].rstrip, #	(-NH3)
				})
			else
				row['theospec_comments'] << "#{line}\n"
			end
		end

	end

	File.open(marshaled_formatted_moda_output, 'wb') {|f| f.write(Marshal.dump(formatted_moda_output))}
end


initial_formatted_moda_output_headers = formatted_moda_output.headers.dup
initial_formatted_moda_output_headers.delete('theospec_comments')
initial_formatted_moda_output_headers.delete('theospec')



##################################################

#	Match theoretical spectra peaks to observed spectra peaks and output to csv/tsv


csvout = CSV.open( "#{moda_base}.#{options[:suffix]}.ion_statistics.tsv",'w',{ col_sep: "\t" })

csvout.puts %w{scan_number rank peptide number_observed_mz number_theoretical_mz max_intensity nb ny mz_b_i mz_y_i theo_mz_b_i theo_mz_y_i b_i intensity_b_i y_i intensity_y_i b_error y_error hyperscore}

hyperout = CSV.open( "#{moda_base}.#{options[:suffix]}.hyper.tsv",'w',{ col_sep: "\t" })
hyperout.puts initial_formatted_moda_output_headers + %w{hyperscore}



#	I'm kinda surprised that I can just read this file again without some type of reset?
formatted_moda_output.each do |row|
	#<CSV::Row "output_index":"2" "spec_index":"19" "observed_MW":"3094.4649" "charge_state":"5" "scan_number":"1535" "rank":"1" "calculated_MW":"3094.4650" "delta_mass":"-0.0001" "score":"46" "probability":"0.3989" "peptide":"R.KTNDKDEKKEDGKQAENDSSNDDKTKK.S" "protein":"sp|Q9BXP5" "pept_position":"301~327" "mod1":"NA" "mod2":"NA" "mod3":"NA" "PlainPeptide":"KTNDKDEKKEDGKQAENDSSNDDKTKK">

	puts "#{row['scan_number']} : #{row['rank']}"
	observed = observed_spectra.find{|s| s['scans'] == row['scan_number'].to_i }		#	FIRST
	number_observed_mz = observed['mzis'].length
	puts "The total number of observed m/z values : #{number_observed_mz}"
	number_theoretical_mz = row['theospec'].length
	puts "The total number of theoretical m/z values : #{number_theoretical_mz}"
	#	observed['mzis'] is an array of pairs. [0]=peak(mz), [1]=intensity(i)
	max_intensity=observed['mzis'].collect{|o|o[1]}.max
	puts "The largest intensity in spectrum whether it is matched or not : #{max_intensity}"



	#	irb(main):001:0> tol=0.5
	#	=> 0.5
	#	irb(main):005:0> (10..0).to_a
	#	=> []	(ranges MUST increase apparently, so ...)
	#	(0..10).to_a.reverse
	#	=> [10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0]
	#	irb(main):011:0> (0..10).to_a.reverse.collect{|v|tol/(2**v)}
	#	=> [0.00048828125, 0.0009765625, 0.001953125, 0.00390625, 0.0078125, 0.015625, 0.03125, 0.0625, 0.125, 0.25, 0.5]

	#	omt (Observed Matched to Theoretical)

	puts "Matching ..."
	omt = []
	mzis = observed['mzis'].dup	#	likely unnecessary as not used after this
#	puts mzis.inspect
	(0..10).to_a.reverse.collect{|v|default_tol/(2**v)}.each do |actual_tol|
#	[default_tol].each do |actual_tol|
		puts "... with tolerance of #{actual_tol}"
		mzis.delete_if do |o|

			indices = row['theospec'].each_index.select{|i|
				( row['theospec'][i][:mz] > ( o[0] - actual_tol ) ) &&
				( row['theospec'][i][:mz] < ( o[0] + actual_tol ) ) }

			if indices.length == 0
				matched = {}
			elsif indices.length == 1
				matched = row['theospec'].delete_at(indices[0])
			else

				puts
				puts "MULTIPLE MATCHES: #{indices.length}"
				puts "MULTIPLE : #{row['theospec'].values_at(*indices)}"
				puts

				# If there are more than one matched theoretical m/z per one observed m/z,
				#	then pick the closest theoretical m/z value to the observed one.
				#	This observed peak is matched.
				#	Ties
				#	Use lowest charge.
				#	Then b ion over y ion
				#	Then higher number ion over lower ion

#MULTIPLE : [
#		{:mz=>636.29876704, :charge_state=>"+1", :ion=>"y6", :d=>"0", :e=>"()"},
#		{:mz=>636.299101391, :charge_state=>"+2", :ion=>"b11", :d=>"0", :e=>"()"}]

#	[{:mz=>230.149917959, :b=>"+1", :ion=>"b2", :d=>"0", :e=>"()", :index=>43, :diff=>0.001203858999986096},
#  {:mz=>230.149917959, :b=>"+1", :ion=>"y2", :d=>"0", :e=>"(-H2O)", :index=>44, :diff=>0.001203858999986096}]

#	Prefer b ion (starting with b) than y ion, so ion => “b7”
#		([{:mz=>831.347173796, :b=>"+1", :ion=>"b7", :d=>"0", :e=>"(-NH3)"}
#	If both are the same type of ions (both b ions and both y ions),
#		then prefer higher number. (For b6 and b7, prefer b7).


				#	by default sorting is ascending... [1,2,3], [a,b,c], ...
				#	Add - before integer for force descending.
				matches = indices.collect{ |i|
					v = row['theospec'][i].dup
					v[:index] = i
					v[:diff] = v[:mz] - o[0]
					v
				}.sort_by{|x| [
					x[:diff].abs,
					x[:charge_state].to_i,
					x[:ion][/([A-z]+)/],
					-x[:ion][/([0-9]+)/].to_i
				] }

				puts "Matches sorted by abs diff, charge state, ion letter, ion number(desc). Selecting the first."
				puts matches.inspect
				puts

				matched = row['theospec'].delete_at(matches.first[:index])

				puts matched.inspect

			end



			omt.delete_if{|x| x[:mz] == o[0] && x[:int] == o[1] }
			omt.push({ mz: o[0], int: o[1], matched: matched })

			!matched.empty?
		end
		puts "Observed match count: #{omt.select{|x|!x[:matched].empty?}.length}"
		puts "Remaining unmatched observed count: #{mzis.length}"
	end
	omt.sort_by!{|o|o[:mz]}

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
	puts "The list of matched observed y-ion intensities divided by the largest intensity in spectrum, and times 100. Let’s call it “Intensity_y_i” for each ith observed peak. – separated by semicolon : #{intensity_y_i}"



#	m/z error = observed m/z minus theoretical m/z

	b_error = bions.collect{|i| i[:mz] - i[:matched][:mz] }.join(';')
	puts "m/z errors between observed and theoretical m/z values for b-ion matches – separated by semicolon : #{b_error}"

	y_error = yions.collect{|i| i[:mz] - i[:matched][:mz] }.join(';')
	puts "m/z errors between observed and theoretical m/z values for y-ion matches – separated by semicolon : #{y_error}"


	mz_b_i = bions.collect{|i| i[:mz] }.join(';')
	puts "Observed matched m/z values (separated by semicolon) for b-ions (col name: mz_b_i) : #{mz_b_i}"


	mz_y_i = yions.collect{|i| i[:mz] }.join(';')
	puts "Observed matched m/z values for y-ions (col name: mz_y_i) : #{mz_y_i}"


	theo_mz_b_i = bions.collect{|i| i[:matched][:mz] }.join(';')
	puts "Theoretical matched m/z values (separated by semicolon) for b-ions (col name: theo_mz_b_i) : #{theo_mz_b_i}"


	theo_mz_y_i = yions.collect{|i| i[:matched][:mz] }.join(';')
	puts "Theoretical matched m/z values for y-ions (col name: theo_mz_y_i) : #{theo_mz_y_i}"







	nbfactorial = (1..nb).inject(1,:*)
	nyfactorial = (1..ny).inject(1,:*)
	b_i_sum = bions.collect{|i|i[:int]}.inject(0){|sum,x| sum + x }
	y_i_sum = yions.collect{|i|i[:int]}.inject(0){|sum,x| sum + x }
	intensity_b_i_sum = bions.collect{|i|i[:int]*100/max_intensity}.inject(0){|sum,x| sum + x }
	intensity_y_i_sum = yions.collect{|i|i[:int]*100/max_intensity}.inject(0){|sum,x| sum + x }
#	hyperscore = Math.log( nbfactorial * nyfactorial * b_i_sum * y_i_sum )
	intensity_hyperscore = Math.log( nbfactorial * nyfactorial * intensity_b_i_sum * intensity_y_i_sum )


	hyperout.puts initial_formatted_moda_output_headers.collect{|x|row[x]} + [ intensity_hyperscore ]

	csvout.puts [row['scan_number'], row['rank'], row['peptide'],
		number_observed_mz, number_theoretical_mz, max_intensity, nb, ny,
		mz_b_i, mz_y_i, theo_mz_b_i, theo_mz_y_i,
		b_i, intensity_b_i, y_i, intensity_y_i, b_error, y_error, intensity_hyperscore]

	puts
	puts
	puts
	puts
end

csvout.close

