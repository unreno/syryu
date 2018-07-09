#!/usr/bin/env ruby

require 'fileutils'

#	running on D: as has more space

source  = 'D:\massive.ucsd.edu\MSV000079053'

target = 'D:\out'

FileUtils.mkdir_p "#{target}" unless File.directory? "#{target}"


#msconvert = 'C:\Program Files\ProteoWizard\ProteoWizard 3.0.18187.b51377ef8\msconvert.exe';
msconvert = 'C:\Program Files\ProteoWizard\ProteoWizard 3.0.18187.b51377ef8\msconvert.exe';

exit unless ARGV.length == 1
bacterium = ARGV[0]



puts bacterium
FileUtils.mkdir_p "#{target}\\#{bacterium}" unless File.directory? "#{target}\\#{bacterium}"

Dir::chdir( "#{source}\\raw\\#{bacterium}" )

Dir["*"].each do |raw|
	puts "-#{raw}"

	puts "-Running msconvert on #{raw}"

	#	lots of quotes are needed
	puts "\"#{msconvert}\" #{raw} --mgf --filter \"msLevel 2\" --filter \"zeroSample removeExtra\" --outdir \"#{target}\\#{bacterium}\"";
	puts `"#{msconvert}" #{raw} --mgf --filter "msLevel 2" --filter "zeroSample removeExtra" --outdir "#{target}\\#{bacterium}"`;

	raw_base = raw.gsub(/\.RAW$/,"")
	mgf = "#{raw_base}.mgf")

	puts "#{source}\\raw\\#{bacterium}\\#{raw}"
	File.delete "#{source}\\raw\\#{bacterium}\\#{raw}"

	Dir::chdir( "#{source}\\sequence\\#{bacterium}" )
	Dir["*"].each do |fasta|

		fasta_base = fasta.gsub(/\.fasta$/,'')
		outdir = "#{target}\\#{bacterium}\\#{fasta_base}"
		FileUtils.mkdir_p "#{outdir}" unless File.directory? "#{outdir}"
		Dir::chdir "#{outdir}"

		puts "--Creating config"
		config = "#{raw_base}.config")
		File.open(config,"w") do |f|

			f.puts "Spectra=#{target}\\#{bacterium}\\#{mgf}"
			f.puts "Fasta=#{source}\\sequence\\#{bacterium}\\#{fasta}"

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

		puts "--Running MODa on #{raw} and #{fasta}"

		out = "#{outdir}/#{raw_base}.out"
		puts "Running java -Xmx10G -jar C:\\ryulab\\moda_v1.51\\moda_v151.jar -i \"#{config}\" -o \"#{out}\""
		puts `java -Xmx10G -jar C:\\ryulab\\moda_v1.51\\moda_v151.jar -i "#{config}" -o "#{out}"`

		puts "gzip #{out}"
		`gzip #{out}`

		puts "#{source}/sequence/#{bacterium}/#{fasta}"
#		File.delete("#{source}/sequence/#{bacterium}/#{fasta}")

	end

	puts "gzip #{target}\\#{bacterium}\\#{mgf}"
	`gzip #{target}\\#{bacterium}\\#{mgf}`

end

