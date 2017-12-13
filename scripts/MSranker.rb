#!/usr/bin/env ruby

require 'csv'
require 'optparse'


def usage(options={})
	puts
	puts "Usage"
	puts
	puts "Uses a date suffix for all created files."
	puts "If going to run multiple times, for whatever reason, you can specify the same suffix with suffix option."
	puts "This way, most of the parsing and the theospec calls won't need rerun."
	puts
	puts "#{$0} [options] <formatted moda output file> <observed spectra mgf file>"
	puts
	puts "Options:"
	puts "	--suffix STRING ...... output file and directory suffixes"
	puts "	--log_count INTEGER .. number of records logged before silencing output"
	puts
	puts "Defaults/current values:"
	puts "	--suffix #{options[:suffix]}"
	puts
	puts "#{$0} out_blind_140521_EOC_MCis_T2_3.txt.OUTPUT.tsv 140521_EOC_MCis_T2_3.SCANS.mgf"
	puts
	exit
end






class Peptide
	attr_reader :peptide
	def initialize(peptide)
		puts "New peptide :#{peptide}:"
		@peptide = peptide
		parts = peptide.split(/\./); #	-.QWLPKD-123EESFLQFKYQALQVP+117.P	(theoretically could include decimals)
		@pre=parts[0]               #	-
		middle=parts[1..-2].join('.') #	QWLPKD-123EESFLQFKYQALQVP+117
		@post=parts[-1]              #	P

		#	parse middle extracting numbers and replacing them with codes passed to theospec
		#@numbers=middle.split(/[A-Z]+/)    #	["",-123,+117]
		@numbers=middle.split(/[A-Z]+/).delete_if{|n|n.empty?}    #	[-123,+117]
		@letters=middle.split(/[0-9.+-]+/) #	[ "QWLPKD","EESFLQFKYQALQVP"]
	end
	def args_for_theospec
		args = ""
#		@numbers.delete_if{|n|n.empty?}.each_with_index do |number,i|
		@numbers.each_with_index do |number,i|
			args << "-c#{i}=#{number} "
		end
		#	-c0=-123 -c1=+117
		args << "#{@pre}."
		@letters.each_index do |i|
			args << @letters[i]
			args << "#{i}" unless( @numbers[i].to_s.empty? )
		end
		#	"QWLPKD1EESFLQFKYQALQVP2"
		args << ".#{@post}"

		puts "Args for theospec :#{args}:"

		#	-c0=-123 -c1=+117 QWLPKD0EESFLQFKYQALQVP1
		args
	end
	def sum_shifts
		#
		#	Instead, can you make it like this?
		#
		#	Please use sum instead of individual modification for any number of modifications.
		#
		#	For example, if there is P+1EP+1TID+1E, then try +3 modification for each amino acid (P+3EPTIDE, PE+3PTIDE, PEP+3TIDE, PEPT+3IDE, so on).
		#
		#	For example, if there is P+1EPTID+1E, then try +2 modification for each amino acid (P+2EPTIDE, PE+2PTIDE, PEP+2TIDE, PEPT+2IDE, so on).
		#
		if @numbers.empty?
			[]
		else
			sum = @numbers.collect{|n|n.to_i}.inject(0,:+)
			out = @letters.join
			@letters.join.length.times.collect do |i|
				"#{@pre}.#{out.dup.insert(i+1,sum.to_s)}.#{@post}"
			end
		end
	end
	def shifts
		#	@letters = ["QEWR","ASDFASDFASDF"]
		#	@numbers = ["-123", "117"]
		#(0...@letters.length).to_a.repeated_permutation(@numbers.length).to_a
		#	=> [[0, 0, 0], [0, 0, 1], [0, 0, 2], [0, 0, 3], [0, 1, 0], [0, 1, 1], [0, 1, 2], [0, 1, 3], ....
		#			[3, 2, 0], [3, 2, 1], [3, 2, 2], [3, 2, 3], [3, 3, 0], [3, 3, 1], [3, 3, 2], [3, 3, 3]]
		(0...@letters.join.length).to_a.repeated_permutation(@numbers.length).collect do |permutation|
#puts permutation.inspect
			out=@letters.join().split('').join(' ').split('').collect{|n| (n==' ')?0:n}.push 0
			#	don't forget to add a trailing 0
			#	=> ["Q", 0, "W", 0, "L", 0, "P", 0, "K", 0, "D", 0, "E", 0, "E", 0, "S", 0, "F", 0, ...
#puts out.inspect
			permutation.each_with_index do |n,i|
#puts [n,i].inspect
				out[2*n+1] += @numbers[i].to_i
			end
#	I was considering putting the "+" back, but it doesn't seem needed.
			"#{@pre}.#{out.delete_if{|n| n == 0}.join()}.#{@post}"
		end
	end
end



class Theospec

	@@theodir = "theospec_output"
	attr_reader :theospecs

	def self.dir=(value)
		@@theodir = value
		Dir.mkdir(@@theodir) unless Dir.exists?(@@theodir)
	end

	def initialize(peptide,charge_state)

		p = Peptide.new(peptide)

#		#	as peptide could include decimals,
#		#	split on decimals to pick off the first and last and then parse the middle
#		parts = peptide.split(/\./); #	-.QWLPKD-123EESFLQFKYQALQVP+117.P	(theoretically could include decimals)
#		a=parts[0]               #	-
#		b=parts[1..-2].join('.') #	QWLPKD-123EESFLQFKYQALQVP+117
#		c=parts[-1]              #	P
#
#		#	parse middle extracting numbers and replacing them with codes passed to theospec
#		numbers=b.split(/[A-Z]+/)    #	["",-123,+117]
#		letters=b.split(/[0-9.+-]+/) #	[ "QWLPKD","EESFLQFKYQALQVP"]
#
#		outpeptide=""
#
#		letters.each_index do |i|
#			outpeptide="#{outpeptide}#{letters[i]}"
#			outpeptide="#{outpeptide}#{i}" unless( numbers[i+1].to_s.empty? )
#		end
#
#		#	"QWLPKD1EESFLQFKYQALQVP2"	(this is actually 0 and 1, not 1 and 2
#
		z=( charge_state.to_i > 5 ) ? 5 : charge_state

		theospec_command="theospec -z#{z}"
#		numbers.delete_if{|n|n.empty?}.each_with_index do |number,i|
#			theospec_command="#{theospec_command} -c#{i}=#{number}"
#		end
#		#	-c0=-123 -c1=+117
#
#		#theospec_output="#{theodir}#{row['scan_number']}.#{row['rank']}.theospec.txt"
		theospec_output="#{@@theodir}#{peptide}.#{charge_state}.theospec.txt"

#		theospec_command="#{theospec_command} #{a}.#{outpeptide}.#{c} > #{theospec_output}"
		theospec_command="#{theospec_command} #{p.args_for_theospec} > #{theospec_output}"
		system(theospec_command) unless File.exists?( theospec_output )

#		row['theospec_comments']=""
#		row['theospec']=[]
		@theospecs=[]
		File.open( theospec_output, 'r' ).each do |line|
			next if line.empty?
			if line =~ /^[[:digit:]]/
				#	Not sure what these all are, nevertheless, split on the semicolon
				##############################################
				#Form of one element:
				#m;z;n;y;c with
				#m = m/z (mass including protons over charge) in concise form
				#z = charge
				#n = name, e.g. `b5'
				#y = currently 0
				#c = neutral changes, e.g. (-H2O)
				##############################################
				#158.092403080(10);+1;y1;0;(-NH3)
				parts=line.split(';')
#				row['theospec'].push({
				@theospecs.push({
					mz: parts[0].to_f,  #	158.092403080(10) to_f -> 158.092403080
					charge_state: parts[1],        #	+1
					ion: parts[2],      #	y1
					d: parts[3],        #	0
					e: parts[4].rstrip, #	(-NH3)
				})
			else
				#row['theospec_comments'] << "#{line}\n"
			end
		end

	end
end	#	Theospec


class PeakMatch

	@@matching_default_tol = 0.5
	attr_reader :observed_matched_to_theoretical

	def self.matching_default_tol=(value)
		@@matching_default_tol = value
	end

	def initialize(working_observed_mzis,working_theospec_mzis)
		#	Attempt to match each observed peak to a theoretical peak

		puts "Matching ..."
		@observed_matched_to_theoretical = []

		#	We will be deleting items so make a duplicate and use that instead.
#	duping unnecessary inside a method or class due to scoping.	YES, I THINK IT STILL IS!
#		working_observed_mzis = observed['mzis'].dup	#	likely unnecessary as not used after this
#		working_theospec_mzis = row['theospec'].dup	#	needed when computing xcorr

		#	puts mzis.inspect

		#	> (10..0).to_a
		#	=> []	(ranges MUST increase apparently, so ...)
		#	> (0..10).to_a.reverse
		#	=> [10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0]
		(0..10).to_a.reverse.collect{|v|
			#	> tol=0.5
			#	(0..10).to_a.reverse.collect{|v|tol/(2**v)}
			#	=> [0.00048828125, 0.0009765625, 0.001953125, 0.00390625, 0.0078125,
			#		0.015625, 0.03125, 0.0625, 0.125, 0.25, 0.5]
			@@matching_default_tol/(2**v)
		}.each do |actual_tol|

			puts "... with tolerance of #{actual_tol}"
			working_observed_mzis.delete_if do |observed_mzi|

				#	Need to do it this way so have indices which can delete
				indices = working_theospec_mzis.each_index.select{|i|
					( working_theospec_mzis[i][:mz] > ( observed_mzi[0] - actual_tol ) ) &&
					( working_theospec_mzis[i][:mz] < ( observed_mzi[0] + actual_tol ) ) }

				if indices.length == 0
					matched = {}
				elsif indices.length == 1
					matched = working_theospec_mzis.delete_at(indices[0])
				else

					puts
					puts "MULTIPLE MATCH COUNT: #{indices.length}"
					puts "MULTIPLE MATCHES : #{working_theospec_mzis.values_at(*indices)}"
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

					#	[{:mz=>230.149917959, :b=>"+1", :ion=>"b2", :d=>"0", :e=>"()", :index=>43, :diff=>0.0012038589},
					#  {:mz=>230.149917959, :b=>"+1", :ion=>"y2", :d=>"0", :e=>"(-H2O)", :index=>44, :diff=>0.0012038589}]

					#	Prefer b ion (starting with b) than y ion, so ion => “b7”
					#		([{:mz=>831.347173796, :b=>"+1", :ion=>"b7", :d=>"0", :e=>"(-NH3)"}
					#	If both are the same type of ions (both b ions and both y ions),
					#		then prefer higher number. (For b6 and b7, prefer b7).


					#	by default sorting is ascending... [1,2,3], [a,b,c], ...
					#	Add - before integer for force descending.
					matches = indices.collect{ |i|
						v = working_theospec_mzis[i].dup
						v[:index] = i		#	NEED this so know what to delete if matches
						v[:diff] = v[:mz] - observed_mzi[0]
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

					matched = working_theospec_mzis.delete_at(matches.first[:index])

					puts matched.inspect

				end

				#	Remove existing non-match if exists.
				@observed_matched_to_theoretical.delete_if do |observed_peak|
					observed_peak[:mz] == observed_mzi[0] &&
																observed_peak[:int] == observed_mzi[1] &&
																observed_peak[:matched] == {}
				end

				@observed_matched_to_theoretical.push({
					mz: observed_mzi[0], int: observed_mzi[1], matched: matched })

				#	Flags the deletion if match is found
				!matched.empty?
			end
			puts "Observed match count: #{@observed_matched_to_theoretical.select{|x|!x[:matched].empty?}.length}"
			puts "Remaining unmatched observed count: #{working_observed_mzis.length}"
		end
		@observed_matched_to_theoretical.sort_by!{|o|o[:mz]}

		puts @observed_matched_to_theoretical.inspect
	end
end	#	PeakMatch


class Hyperscore

	attr_reader :nb, :ny, :b_i, :y_i, :intensity_b_i, :intensity_y_i, :b_error, :y_error,
		:mz_b_i, :mz_y_i, :theo_mz_b_i, :theo_mz_y_i, :hyperscore

	def initialize(bions,yions,max_intensity)
		puts "bions :#{bions}:"
		puts "yions :#{yions}:"
		puts "max_intensity : #{max_intensity}"

		@nb = bions.length
		puts "# of observed b-ion m/z values matched (b-ion is starting with b in theoretical spec; for example, b2, b3, b4, ….) Let’s call it as “Nb” : #{nb}"

		@ny = yions.length
		puts "# of observed y-ion m/z values matched (y-ion is starting with y in theoretical spec; for example, y2, y5,…) Let’s call it as “Ny” : #{ny}"

		@b_i = bions.collect{|i|i[:int]}.join(';')
		puts "The list of matched observed b-ion intensities in the spectrum (in mgf files). – separated by semicolon : #{b_i}"

		@intensity_b_i = bions.collect{|i|i[:int]*100/max_intensity}.join(';')
		puts "The list of matched observed b-ion intensities divided by the largest intensity in spectrum, and times 100. Let’s call it “Intensity_b_i” for each ith observed peak.  – separated by semicolon : #{intensity_b_i}"

		@y_i = yions.collect{|i|i[:int]}.join(';')
		puts "The list of matched observed y-ion intensities (in mgf files) – separated by semicolon : #{y_i}"

		@intensity_y_i = yions.collect{|i|i[:int]*100/max_intensity}.join(';')
		puts "The list of matched observed y-ion intensities divided by the largest intensity in spectrum, and times 100. Let’s call it “Intensity_y_i” for each ith observed peak. – separated by semicolon : #{intensity_y_i}"



		#	m/z error = observed m/z minus theoretical m/z

		@b_error = bions.collect{|i| i[:mz] - i[:matched][:mz] }.join(';')
		puts "m/z errors between observed and theoretical m/z values for b-ion matches – separated by semicolon : #{b_error}"

		@y_error = yions.collect{|i| i[:mz] - i[:matched][:mz] }.join(';')
		puts "m/z errors between observed and theoretical m/z values for y-ion matches – separated by semicolon : #{y_error}"

		@mz_b_i = bions.collect{|i| i[:mz] }.join(';')
		puts "Observed matched m/z values (separated by semicolon) for b-ions (col name: mz_b_i) : #{mz_b_i}"

		@mz_y_i = yions.collect{|i| i[:mz] }.join(';')
		puts "Observed matched m/z values for y-ions (col name: mz_y_i) : #{mz_y_i}"

		@theo_mz_b_i = bions.collect{|i| i[:matched][:mz] }.join(';')
		puts "Theoretical matched m/z values (separated by semicolon) for b-ions (col name: theo_mz_b_i) : #{theo_mz_b_i}"

		@theo_mz_y_i = yions.collect{|i| i[:matched][:mz] }.join(';')
		puts "Theoretical matched m/z values for y-ions (col name: theo_mz_y_i) : #{theo_mz_y_i}"

		#	Set initial value to 1 just in case N is 0, otherwise returns nil
		#	> (1..0).inject(:*)
		#	=> nil
		#	> (1..0).inject(1,:*)
		#	=> 1

		nbfactorial = (1..nb).inject(1,:*)
		puts "Nb! : #{nbfactorial}"

		nyfactorial = (1..ny).inject(1,:*)
		puts "Ny! : #{nyfactorial}"

		intensity_b_i_sum = bions.collect{|i|i[:int]*100/max_intensity}.inject(0,:+)	#	sum
		puts "intensity_b_i_sum : #{intensity_b_i_sum}"

		intensity_y_i_sum = yions.collect{|i|i[:int]*100/max_intensity}.inject(0,:+)	#	sum
		puts "intensity_y_i_sum : #{intensity_y_i_sum}"

		pre_hyperscore_product = nbfactorial * nyfactorial * intensity_b_i_sum * intensity_y_i_sum

		@hyperscore = if pre_hyperscore_product == 0
			puts "Pre Hyperscore Product is 0 so Hyperscore will be -Infinity. Setting it to -9999 instead."
			-9999
		else
			Math.log( pre_hyperscore_product )
		end

		puts "Hyperscore : #{hyperscore}"
	end

end	#	Hyperscore



class Xcorr
	attr_reader :xcorr
	def initialize(x0,y0,xcorr_max_mass,xcorr_bin_size)
		#	The following is how to calculate the second score “xcorr” for MODa.
		#	1. For each sequence in MODa, compute theoretical m/z (you already have codes)
		#	2. For each input in #1, you can find the corresponding scan in mgf file. (you already have codes)

		#	3. Binning theoretical m/z from #1 using bin.size=2 (by default) and max.mass=2500 (by default)
		#	a. Create a vector (or array) x0 with its size max.mass/bin.size=2500/2=1250
		#	b. A vector will contain zero or one. If there exists at least one peak corresponding to bin location,
		#		then x0[bin.location]=1. Otherwise, x0[bin.location]=0.
		#		(Not sure if a programming language you use has array (or vector) number starting from 0 or 1,
		#		but assuming it starts from 1.)
		#	Note that x0[1] represents peaks with their 0 < m/z value <= bin.size.
		#	Note that x0[2] represents peaks with their bin.size < m/z value <= 2*bin.size.
		#	Thus, if we have only theoretical m/z value 2.5 (unrealistic example), then all are zero except x0[2]=1.
		#	If we have only theoretical m/z values 2.5 and 2.3, then we have the same vector with all zero except x0[2]=1.
		#	If there is peak larger than max.mass, ignore this (just give a warning to increase max.mass).

		puts "x0 :#{x0}:"

		#	4. Binning observed m/z from #2 using the same bin.size and max.mass as above.
		#	a. First, observed intensity (named as intensity.one) will range from 0 to 1. So for hyperscore,
		#		we divided peak intensity by the largest peak times 100. But for Xcorr, we will not multiply by
		#		100. Thus, peak intensity/the largest peak intensity is the intensity we will use. Let’s call
		#		this intensity as intensity.one.

		#	b. Create a vector (or array) (same as 2a) y0 with the same bin.size and max.mass.
		#	c. A vector (or array) will contain total intensities of peaks belong to that bin.
		#	For example y0[1] is the total intensities with their peak m/z between 0 and bin.size (including bin.size).
		#	If we have only observed peaks with their m/z values 2.5 (with intensity 0.5) and m/z value 2.3
		#		(with intensity 0.4), then y0[2]=0.9 and others are zero.

		puts "y0 :#{y0}:"

		#	5. Finally, calculate Xcorr. Complicated! But let me explain.
		#	a. sum.off and y.dash are vectors (or arrays) with the same size as ones created in #2 and #3. Create these.
		#	b. Initially, all zero in sum.off vector and y.dash.

		xcorr_sum_off = Array.new(xcorr_max_mass/xcorr_bin_size,0)
		xcorr_y_dash  = Array.new(xcorr_max_mass/xcorr_bin_size,0)

		#	c. Calculate sum.off and y.dash. I assume that an array (or vector) starts from location one, again.
		#	For j in 1:(size of sum.off){
		#		For k in -75:75 {
		#			If (k is not 0){
		#				If (j+k >0 AND j+k <= (size of y0)){
		#					sum.off[j] = sum.off[j] + y0[j+k]
		#				}
		#			}
		#		}
		#		y.dash[j] = y0[j] – (sum.off[j] / 150)
		#	}

		(0...(xcorr_max_mass/xcorr_bin_size)).to_a.each do |j| # 0-1249
			(-75..75).to_a.each do |k|	#	-75-75
				if k != 0 then
					#	make sure that the subscript will be within the bounds of these arrays.
					if( j+k >= 0 && j+k < y0.length ) then
						xcorr_sum_off[j] = xcorr_sum_off[j] + y0[j+k]
					end
				end
			end
			xcorr_y_dash[j] = y0[j] - ( xcorr_sum_off[j] / 150 )
		end

		puts "xcorr_sum_off: #{xcorr_sum_off}"
		puts "xcorr_y_dash: #{xcorr_y_dash}"



		#	d. Compute dot product between x0 and y.dash. If your program language does not have dot
		#		product function, then you can do the following.
		#	Initialize xcorr=0
		#	For j in (1: size of x0){
		#		xcorr=xcorr+ x0[j] * y.dash[j]
		#	}
		#	Here, * is a regular multiplication.

		@xcorr = 0

		#	could do the following, but oddly benchmarks show it slower.
		#	require 'matrix'
		#	Vector[*x0].dot(Vector[*xcorr_y_dash])

		(0...(xcorr_max_mass/xcorr_bin_size)).to_a.each do |j| # 0-1249
			@xcorr += ( x0[j] * xcorr_y_dash[j] )
		end

		puts "xcorr :#{@xcorr}"
	end

end	#	Xcorr



######################################################################



stdout_redirected = false
original_stdout = $stdout.clone
#	$stdout.reopen("/dev/null", "a")


#	Must be called before option parsing as they remove the items.
puts "Command: #{$0} #{$*.join(' ')}"


#options = { suffix: 'TESTING'}
options = { suffix: Time.now.strftime("%Y%m%d%H%M%S") }
OptionParser.new do |opt|
	opt.on('--suffix TESTING') { |o| options[:suffix] = o }
	opt.on('--log_count 10') { |o| options[:log_count] = o }
end.parse!

#puts options

puts "Using suffix '#{options[:suffix]}'"

usage(options) if ARGV.length < 2


#	This gets ugly if the filename begins with a dot and has no extension.
#
#in_moda_filename_without_extension = File.basename( ARGV[0], ".*" )	#	loses path
#in_moda_filename_without_extension = ARGV[0].sub(/#{File.extname(ARGV[0])}$/,'')
in_moda_filename_without_extension = ARGV[0].chomp(File.extname(ARGV[0]))

#	Never actually used, but just in case
#in_mgf_filename_without_extension = File.basename( ARGV[1], ".*" )	#	loses path
#in_mgf_filename_without_extension = ARGV[1].sub(/#{File.extname(ARGV[1])}$/,'')
#in_mgf_filename_without_extension = ARGV[1].chomp(File.extname(ARGV[1]))


#theodir="#{ARGV[0]}_theospec_#{options[:suffix]}/"
#Dir.mkdir(theodir) unless Dir.exists?(theodir)
Theospec.dir="#{ARGV[0]}_theospec_#{options[:suffix]}/"


#	Let’s set the default m/z error as tol=0.5. (this is used when matching observed to theoretical)
#PeakMatch.matching_default_tol = 0.5

#	Default settings when computing xcorr
xcorr_bin_size = 1	#	2
xcorr_max_mass = 2500

##################################################

#	Read, parse and store observed spectra.

#	Sample snippet.
#	BEGIN IONS
#	TITLE=140521_EOC_MCis_T2_3.11.11.2 File:"140521_EOC_MCis_T2_3.raw", NativeID:"controllerType=0 controllerNumber=1 scan=11"
#	RTINSECONDS=3.04065084
#	PEPMASS=573.833417843252 7326.4638671875
#	SCANS=11
#	CHARGE=2+
#	138.3958823 599.2578735352
#	149.0234809 10632.353515625
#	159.8794218 661.0111694336
#	205.0861157 2093.7409667969
#	251.8871732 747.2424926758
#	252.8723301 883.7633666992
#	273.4169172 585.7670898438
#	302.3671104 577.3031005859
#	307.9346258 593.6673583984
#	660.8964981 659.2720336914
#	670.1807365 640.5803833008
#	677.2812549 711.8297729492
#	END IONS

#	After running the data is saved as a marshal which can be
#	reused, if desired, by specifying the same suffix.
#	This saves time, particularly useful during development.

marshaled_observed_spectra_filename="#{ARGV[1]}.#{options[:suffix]}.marshal"
if File.exists?( marshaled_observed_spectra_filename )
	observed_spectra=Marshal.load(File.binread(marshaled_observed_spectra_filename))
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
			puts "Processing MGF block for scans #{$1.to_i}"
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
	File.open(marshaled_observed_spectra_filename, 'wb') {|f| f.write(Marshal.dump(observed_spectra))}
end




##################################################

#	Read and parse moda output, then compute and store theoretical spectra.

#	Sample "head -2"
#	output_index	spec_index	observed_MW	charge_state	scan_number	rank	calculated_MW	delta_mass	score	probability	peptide	protein	pept_position	mod1	mod2	mod3	PlainPeptide
#	2	19	3094.4649	5	1535	1	3094.4650	-0.0001	46	0.3989	R.KTNDKDEKKEDGKQAENDSSNDDKTKK.S	sp|Q9BXP5	301~327	NA	NA	NA	KTNDKDEKKEDGKQAENDSSNDDKTKK
#	....

#	After running the data is saved as a marshal which can be
#	reused, if desired, by specifying the same suffix.
#	This saves time, particularly useful during development.

marshaled_formatted_moda_filename="#{ARGV[0]}.#{options[:suffix]}.marshal"
if File.exists?( marshaled_formatted_moda_filename )
	formatted_moda_output=Marshal.load(File.binread(marshaled_formatted_moda_filename))
else
	formatted_moda_output = CSV.read( ARGV[0],'rb',{ col_sep: "\t", headers: true })
	formatted_moda_output.each do |row|

		puts "Processing MODa row #{row['scan_number']}:#{row['rank']}"

	#<CSV::Row "output_index":"2" "spec_index":"19" "observed_MW":"3094.4649" "charge_state":"5" "scan_number":"1535" "rank":"1" "calculated_MW":"3094.4650" "delta_mass":"-0.0001" "score":"46" "probability":"0.3989" "peptide":"R.KTNDKDEKKEDGKQAENDSSNDDKTKK.S" "protein":"sp|Q9BXP5" "pept_position":"301~327" "mod1":"NA" "mod2":"NA" "mod3":"NA" "PlainPeptide":"KTNDKDEKKEDGKQAENDSSNDDKTKK">

		row['theospec'] = Theospec.new( row['peptide'], row['charge_state'] ).theospecs

	end

	File.open(marshaled_formatted_moda_filename, 'wb') {|f| f.write(Marshal.dump(formatted_moda_output))}
end


#initial_formatted_moda_output_headers = formatted_moda_output.headers.dup
#initial_formatted_moda_output_headers.delete('theospec_comments')
#initial_formatted_moda_output_headers.delete('theospec')
#	output_index	spec_index	observed_MW	charge_state	scan_number	rank	calculated_MW	delta_mass	score	probability	peptide	protein	pept_position	mod1	mod2	mod3	PlainPeptide

initial_formatted_moda_output_headers = %w{output_index	spec_index	observed_MW	charge_state	scan_number	rank	calculated_MW	delta_mass	score	probability	peptide	protein	pept_position	mod1	mod2	mod3	PlainPeptide}




##################################################

#	Match theoretical spectra peaks to observed spectra peaks and output to csv/tsv

ion_statistics_columns = %w{scan_number rank peptide number_observed_mz number_theoretical_mz max_intensity nb ny mz_b_i mz_y_i theo_mz_b_i theo_mz_y_i b_i intensity_b_i y_i intensity_y_i b_error y_error hyperscore}

ion_statistics_file = CSV.open(
	"#{in_moda_filename_without_extension}.#{options[:suffix]}.ion_statistics.tsv",
	'w',{ col_sep: "\t" })

ion_statistics_file.puts ion_statistics_columns

basic_ion_statistics_file = CSV.open(
	"#{in_moda_filename_without_extension}.#{options[:suffix]}.basic_ion_statistics.tsv",
	'w',{ col_sep: "\t" })

basic_ion_statistics_file.puts ion_statistics_columns

extended_moda_file = CSV.open(
	"#{in_moda_filename_without_extension}.#{options[:suffix]}.hyper.tsv",
	'w',{ col_sep: "\t" })

extended_moda_file.puts initial_formatted_moda_output_headers + %w{hyperscore xcorr hyperscore.basic xcorr.basic peptide.hyperscore.basic.with.shift hyperscore.basic.with.shift peptide.xcorr.basic.with.shift xcorr.basic.with.shift}






#	I'm kinda surprised that I can just read this file again without some type of reset?
#	I guess having been fully read makes it ok.
formatted_moda_output.each_with_index do |row,record_number|

	if !stdout_redirected && options.has_key?(:log_count) && record_number > options[:log_count].to_i
		puts "Record number(#{record_number}) has exceeded requested log count(#{options[:log_count]})."
		puts "Redirecting the rest of the output to /dev/null."
		stdout_redirected = true
		$stdout.reopen("/dev/null", "a")
	end

	#	Not sure if "row" would preserve order so ...
	extended_moda_file_array = initial_formatted_moda_output_headers.collect{|x|row[x]}

	ion_statistics_file_array = [ row['scan_number'], row['rank'], row['peptide'] ]

	#<CSV::Row "output_index":"2" "spec_index":"19" "observed_MW":"3094.4649" "charge_state":"5" "scan_number":"1535" "rank":"1" "calculated_MW":"3094.4650" "delta_mass":"-0.0001" "score":"46" "probability":"0.3989" "peptide":"R.KTNDKDEKKEDGKQAENDSSNDDKTKK.S" "protein":"sp|Q9BXP5" "pept_position":"301~327" "mod1":"NA" "mod2":"NA" "mod3":"NA" "PlainPeptide":"KTNDKDEKKEDGKQAENDSSNDDKTKK">

	puts "SCAN:MODa_RANK => #{row['scan_number']} : #{row['rank']}"
	#	find returns only FIRST, but there should only be 1 match
	observed = observed_spectra.find{|s| s['scans'] == row['scan_number'].to_i }
	ion_statistics_file_array.push number_observed_mz = observed['mzis'].length
	puts "The total number of observed m/z values : #{number_observed_mz}"

	ion_statistics_file_array.push number_theoretical_mz = row['theospec'].length
	puts "The total number of theoretical m/z values : #{number_theoretical_mz}"

	#	observed['mzis'] is an array of pairs. [0]=peak(mz), [1]=intensity(i)
	ion_statistics_file_array.push max_intensity = observed['mzis'].collect{|o| o[1] }.max
	puts "The largest intensity in spectrum whether it is matched or not : #{max_intensity}"



	p=PeakMatch.new( observed['mzis'].dup, row['theospec'].dup )

	puts
	puts "HYPERSCORE CALCULATION (ALL)"

	h1 = Hyperscore.new(
		p.observed_matched_to_theoretical.select{|o| o[:matched][:ion] =~ /^b/ },
		p.observed_matched_to_theoretical.select{|o| o[:matched][:ion] =~ /^y/ },
		max_intensity)

	puts
	puts "HYPERSCORE CALCULATION (BASIC EMPTY PARENS)"

	h2 = Hyperscore.new(
		p.observed_matched_to_theoretical.select{|o| o[:matched][:ion] =~ /^b/ && o[:matched][:e] == "()" },
		p.observed_matched_to_theoretical.select{|o| o[:matched][:ion] =~ /^y/ && o[:matched][:e] == "()" },
		max_intensity)


	puts
	puts "XCORR CALCULATION"

	#
	#	working_theospec_mzis should have less as matches deleted? Hmm. Why am I using? Just logging. Testing.
	#	I ran with both and no difference in log output. Trying with full log. Difference!
	#	If different, row is always bigger.
	#
	#	greater_than_max_mass = working_theospec_mzis.select{|t| t[:mz] > xcorr_max_mass }
	greater_than_max_mass = row['theospec'].select{|t| t[:mz] > xcorr_max_mass }

	if greater_than_max_mass.length > 0 then
		puts "WARNING: #{greater_than_max_mass.length} theoretical peaks greater than #{xcorr_max_mass}"
		puts "greater_than_max_mass:#{greater_than_max_mass}"
	end

	x1 = Xcorr.new(
		(0...(xcorr_max_mass/xcorr_bin_size)).to_a.collect{|bin|	#	0-1249
			( row['theospec'].select{ |theoretical_peak|
					theoretical_peak[:mz] > bin*xcorr_bin_size &&
					theoretical_peak[:mz] <= (bin+1)*xcorr_bin_size
				}.empty? ) ? 0 : 1
		},
		(0...(xcorr_max_mass/xcorr_bin_size)).to_a.collect{|bin|	#	0-1249
			p.observed_matched_to_theoretical.select{|observed_peak|
				observed_peak[:mz] > bin*xcorr_bin_size && observed_peak[:mz] <= (bin+1)*xcorr_bin_size
			}.collect{|observed_peak|
				observed_peak[:int]/max_intensity
			}.inject(0,:+) #	sums array (set 0 as initial value for empty arrays otherwise returns nil)
		},xcorr_max_mass,xcorr_bin_size )

	x2 = Xcorr.new(
		(0...(xcorr_max_mass/xcorr_bin_size)).to_a.collect{|bin|	#	0-1249
			( row['theospec'].select{ |theoretical_peak|
					theoretical_peak[:e] == "()" &&
					theoretical_peak[:mz] > bin*xcorr_bin_size &&
					theoretical_peak[:mz] <= (bin+1)*xcorr_bin_size
				}.empty? ) ? 0 : 1
		},
		(0...(xcorr_max_mass/xcorr_bin_size)).to_a.collect{|bin|	#	0-1249
			p.observed_matched_to_theoretical.select{|observed_peak|
				observed_peak[:mz] > bin*xcorr_bin_size && observed_peak[:mz] <= (bin+1)*xcorr_bin_size
			}.collect{|observed_peak|
				observed_peak[:int]/max_intensity
			}.inject(0,:+) #	sums array (set 0 as initial value for empty arrays otherwise returns nil)
		},xcorr_max_mass,xcorr_bin_size )


	#	1.       score name: hyperscore.basic
	#	This is the same as hyperscore except their theoretical mass list is shorter.
	#	Before, we considered all the theoretical masses. Now for this score, I want to
	#		include any ions with empty parentheses.
	#	For example, we include 98.0600403189(54);+1;b1;0;().
	#	But we don’t include a mass with something in (). For example, (-NH3) or (-H2O).

	#	2.       score name: xcorr.basic
	#	Calculate xcorr with the same mass list used in #1.

	#	3.       Score name: hyperscore.basic.with.shift
	#	Additional column: peptide.hyperscore.basic.with.shift
	#	Top hyperscore of peptides with varying mass shift. We will use mass list only including
	#		empty () like #1 and #2.
	#	For example, let say we have “K.PEP+10TID+85E.R”
	#	We can have   K.P+95EPTIDE.R, K.P+10E+85PTIDE.R, K.P+10EP+85TIDE.R, K.P+10EPT+85IDE.R,
	#		K.P+10EPTI+85DE.R, K.P+10EPTID+85E.R, K.P+10EPTID+85E.R, K.P+10EPTIDE+85.R,
	#	K.P+85E+10PTIDE.R, K.PE+95PTIDE.R, K.PE+10P+85TIDE.R, K.PE+10PT+85IDE.R, so on….
	#	Here, we try all the mass shift using the same letters. Calculate hyperscore for each sequence.
	#		(Keep all hyperscores in processing.txt files)
	#	Then, take the highest hyperscore and record the corresponding peptide sequence with modification.
	#		If there is a tie, then we can record their sequences separated by semicolon.

	#	4.       Score name: xcorr.basic.with.shift
	#	Additional column: peptide.xcorr.basic.with.shift
	#	Same as #3 except that we use xcorr score. We will use mass list only including empty () like #1 and #2.


	hshiftmax = 0;
	hshiftpeptide = "";
	xshiftmax = 0;
	xshiftpeptide = "";

	puts "A Whole lotta Shifting (may be) about to happen ..."

	uniq_shifts = Peptide.new(row['peptide']).shifts.uniq  - [row['peptide']]

	puts "Uniq Shifts :#{uniq_shifts}:"

	uniq_shifts.each do |shifted_peptide|

		puts "#{row['peptide']} shifted to #{shifted_peptide}"

		theospec_shift = Theospec.new( shifted_peptide, row['charge_state'] ).theospecs

		pshift = PeakMatch.new( observed['mzis'].dup, theospec_shift.dup )

		puts
		puts "HYPERSCORE CALCULATION (BASIC EMPTY PARENS) WITH SHIFT"

		hshift = Hyperscore.new(
			pshift.observed_matched_to_theoretical.select{|o| o[:matched][:ion] =~ /^b/ && o[:matched][:e] == "()" },
			pshift.observed_matched_to_theoretical.select{|o| o[:matched][:ion] =~ /^y/ && o[:matched][:e] == "()" },
			max_intensity)

		if hshift.hyperscore == hshiftmax then
			hshiftpeptide += ";#{shifted_peptide}"
		elsif hshift.hyperscore > hshiftmax
			hshiftmax = hshift.hyperscore
			hshiftpeptide = shifted_peptide
		end

		puts
		puts "XCORR CALCULATION (BASIC EMPTY PARENS) WITH SHIFT"

		xshift = Xcorr.new(
			(0...(xcorr_max_mass/xcorr_bin_size)).to_a.collect{|bin|	#	0-1249
				( row['theospec'].select{ |theoretical_peak|
						theoretical_peak[:e] == "()" &&
						theoretical_peak[:mz] > bin*xcorr_bin_size &&
						theoretical_peak[:mz] <= (bin+1)*xcorr_bin_size
					}.empty? ) ? 0 : 1
			},
			(0...(xcorr_max_mass/xcorr_bin_size)).to_a.collect{|bin|	#	0-1249
				pshift.observed_matched_to_theoretical.select{|observed_peak|
					observed_peak[:mz] > bin*xcorr_bin_size && observed_peak[:mz] <= (bin+1)*xcorr_bin_size
				}.collect{|observed_peak|
					observed_peak[:int]/max_intensity
				}.inject(0,:+) #	sums array (set 0 as initial value for empty arrays otherwise returns nil)
			},xcorr_max_mass,xcorr_bin_size )

		if xshift.xcorr == xshiftmax then
			xshiftpeptide += ";#{shifted_peptide}"
		elsif xshift.xcorr > xshiftmax
			xshiftmax = xshift.xcorr
			xshiftpeptide = shifted_peptide
		end

	end #	Peptide.new(row['peptide']).shifts.each do |shifted_peptide|

	extended_moda_file.puts extended_moda_file_array + [
		h1.hyperscore, x1.xcorr, h2.hyperscore, x2.xcorr, hshiftpeptide, hshiftmax, xshiftpeptide, xshiftmax ]

	#	for initial hyperscore
	ion_statistics_file.puts ion_statistics_file_array + [ h1.nb, h1.ny,
		h1.mz_b_i, h1.mz_y_i, h1.theo_mz_b_i, h1.theo_mz_y_i,
		h1.b_i, h1.intensity_b_i, h1.y_i, h1.intensity_y_i, h1.b_error, h1.y_error, h1.hyperscore]

	#	for new basic hyperscore
	basic_ion_statistics_file.puts ion_statistics_file_array + [ h2.nb, h2.ny,
		h2.mz_b_i, h2.mz_y_i, h2.theo_mz_b_i, h2.theo_mz_y_i,
		h2.b_i, h2.intensity_b_i, h2.y_i, h2.intensity_y_i, h2.b_error, h2.y_error, h2.hyperscore]

	puts
	puts
	puts
	puts
end

extended_moda_file.close
ion_statistics_file.close
basic_ion_statistics_file.close

