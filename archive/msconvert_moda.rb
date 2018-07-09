#!/usr/bin/env ruby

require 'fileutils'

source  = "/Users/jake/massive.ucsd.edu/MSV000079053"

out_base = "/Users/jake/out"

FileUtils.mkdir_p "#{out_base}" unless File.directory? "#{out_base}"


#msconvert = 'C:\Program Files\ProteoWizard\ProteoWizard 3.0.18187.b51377ef8\msconvert.exe';
msconvert = '/Program Files/ProteoWizard/ProteoWizard 3.0.18187.b51377ef8/msconvert.exe';

exit unless ARGV.length == 1
bacterium = ARGV[0]



puts bacterium

Dir::chdir( "#{source}/raw/#{bacterium}" )

Dir["*"].each do |raw|
	puts "-#{raw}"

	FileUtils.mkdir_p "#{out_base}/#{bacterium}" unless File.directory? "#{out_base}/#{bacterium}"

	puts "-Running msconvert on #{raw}"

	#	lots of quotes are needed
	puts "\"#{msconvert}\" #{raw} --mgf --filter \"msLevel 2\" --filter \"zeroSample removeExtra\" --outdir \"#{out_base}/#{bacterium}\"";
	puts `"#{msconvert}" #{raw} --mgf --filter "msLevel 2" --filter "zeroSample removeExtra" --outdir "#{out_base}/#{bacterium}"`;

	mgf = raw.gsub(/RAW$/,"mgf")

	puts "#{source}/raw/#{bacterium}/#{raw}"
	File.delete "#{source}/raw/#{bacterium}/#{raw}"

	Dir::chdir( "#{source}/sequence/#{bacterium}" )
	Dir["*"].each do |sequence|

		outdir = "#{out_base}/#{bacterium}/#{sequence.gsub(/\.fasta$/,'')}"
		FileUtils.mkdir_p "#{outdir}" unless File.directory? "#{outdir}"
		Dir::chdir "#{outdir}"

		puts "--Creating config"
		config = mgf.gsub(/mgf$/,"config")
		File.open(config,"w") do |f|

			f.puts "Spectra=#{out_base}/#{bacterium}/#{mgf}"
			f.puts "Fasta=#{source}/sequence/#{bacterium}/#{sequence}"

			f.puts "Instrument=ESI-TRAP"
			f.puts "PeptTolerance=0.5"
			f.puts "AutoPMCorrection=1"
			f.puts "FragTolerance="
			f.puts "BlindMode=2"
			f.puts "MinModSize="
			f.puts "MaxModSize="
			f.puts "Enzyme=Trypsin, KR/C"
			f.puts "MissedCleavage=2"
			f.puts "Protocol=NONE"
			f.puts "HighResolution=OFF"
		end

		puts "--Running MODa on #{raw} and #{sequence}"

		out = "#{outdir}/#{mgf.gsub(/mgf$/,"out")}"
		puts "Running java -Xmx5000M -jar /ryulab/moda_v1.51/moda_v151.jar -i \"#{config}\" -o \"#{out}\""
		puts `java -Xmx5000M -jar /ryulab/moda_v1.51/moda_v151.jar -i "#{config}" -o "#{out}"`

#		puts "#{source}/sequence/#{bacterium}/#{sequence}"
#		File.delete("#{source}/sequence/#{bacterium}/#{sequence}")

	end

	puts "gzip #{source}/raw/#{bacterium}/#{mgf}"
	`gzip --best #{source}/raw/#{bacterium}/#{mgf}`

end

