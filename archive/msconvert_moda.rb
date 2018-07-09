#!/usr/bin/env ruby

require 'fileutils'

#	running on D:

source  = 'D:\massive.ucsd.edu\MSV000079053'

out_base = 'D:\out'

FileUtils.mkdir_p "#{out_base}" unless File.directory? "#{out_base}"


#msconvert = 'C:\Program Files\ProteoWizard\ProteoWizard 3.0.18187.b51377ef8\msconvert.exe';
msconvert = 'C:\Program Files\ProteoWizard\ProteoWizard 3.0.18187.b51377ef8\msconvert.exe';

exit unless ARGV.length == 1
bacterium = ARGV[0]



puts bacterium
out_bacterium = "#{out_base}\\#{bacterium}"
FileUtils.mkdir_p "#{out_bacterium}" unless File.directory? "#{out_bacterium}"


Dir::chdir( "#{source}\\raw\\#{bacterium}" )

Dir["*"].each do |raw|
	puts "-#{raw}"

	puts "-Running msconvert on #{raw}"

	#	lots of quotes are needed
	puts "\"#{msconvert}\" #{raw} --mgf --filter \"msLevel 2\" --filter \"zeroSample removeExtra\" --outdir \"#{out_bacterium}\"";
	puts `"#{msconvert}" #{raw} --mgf --filter "msLevel 2" --filter "zeroSample removeExtra" --outdir "#{out_bacterium}"`;

	mgf = raw.gsub(/RAW$/,"mgf")
	out_mgf = "#{out_bacterium}\\#{mgf}"

	puts "#{source}\\raw\\#{bacterium}\\#{raw}"
#	File.delete "#{source}\\raw\\#{bacterium}\\#{raw}"

	Dir::chdir( "#{source}\\sequence\\#{bacterium}" )
	Dir["*"].each do |sequence|

		outdir = "#{out_bacterium}\\#{sequence.gsub(/\.fasta$/,'')}"
		FileUtils.mkdir_p "#{outdir}" unless File.directory? "#{outdir}"
		Dir::chdir "#{outdir}"

		puts "--Creating config"
		config = mgf.gsub(/mgf$/,"config")
		File.open(config,"w") do |f|

			f.puts "Spectra=#{out_mgf}"
			f.puts "Fasta=#{source}\\sequence\\#{bacterium}\\#{sequence}"

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
		puts "Running java -Xmx10G -jar C:\\ryulab\\moda_v1.51\\moda_v151.jar -i \"#{config}\" -o \"#{out}\""
		puts `java -Xmx10G -jar C:\\ryulab\\moda_v1.51\\moda_v151.jar -i "#{config}" -o "#{out}"`

		puts "gzip #{out}"
		`gzip --best #{out}`

#		puts "#{source}/sequence/#{bacterium}/#{sequence}"
#		File.delete("#{source}/sequence/#{bacterium}/#{sequence}")

	end

	puts "gzip #{out_mgf}"
	`gzip --best #{out_mgf}`

end

