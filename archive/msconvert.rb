#!/usr/bin/env ruby

#	\Ruby25-x64\bin\ruby.exe

require 'parallel'

require 'net/ftp'
#https://ruby-doc.org/stdlib-2.4.0/libdoc/net/ftp/rdoc/Net/FTP.html

require 'fileutils'

ftp = Net::FTP.new('massive.ucsd.edu')
ftp.login
puts "Logged in"

local_root  = "/Users/jakewendt/"
remote_root = "/"
base = "MSV000079053"
local_base  = "#{local_root}#{base}"
remote_base = "#{remote_root}#{base}"

FileUtils.mkdir_p "#{local_base}/raw" unless File.directory? "#{local_base}/raw"
FileUtils.mkdir_p "#{local_base}/sequence" unless File.directory? "#{local_base}/sequence"
FileUtils.mkdir_p "#{local_base}/out" unless File.directory? "#{local_base}/out"


puts "chdir #{remote_base}/raw"
ftp.chdir("#{remote_base}/raw")



#msconvert = 'C:\Program Files\ProteoWizard\ProteoWizard 3.0.18187.b51377ef8\msconvert.exe';
msconvert = '/Program Files/ProteoWizard/ProteoWizard 3.0.18187.b51377ef8/msconvert.exe';

#
##	parallelize like ... syryu/scripts/modification_locator.rb
#



ftp.nlst.each do |bacterium|
#Parallel.each(ftp.nlst) do |bacterium|
	puts bacterium

	ftp.chdir("#{remote_base}/raw/#{bacterium}")

#	Parallel.each(ftp.nlst) do |raw|
	ftp.nlst.each do |raw|
		puts "-#{raw}"

		FileUtils.mkdir_p "#{local_base}/out/#{bacterium}" unless File.directory? "#{local_base}/out/#{bacterium}"

		FileUtils.mkdir_p "#{local_base}/raw/#{bacterium}" unless File.directory? "#{local_base}/raw/#{bacterium}"
		Dir::chdir "#{local_base}/raw/#{bacterium}"

		ftp.chdir("#{remote_base}/raw/#{bacterium}")

		if File.exists? raw
			puts "-#{raw} already exists"
		else
			puts "-Getting #{raw}"
			ftp.get(raw)
		end



		puts "-Running msconvert on #{raw}"

		#	lots of quotes are needed
		puts "\"#{msconvert}\" #{raw} --mgf --filter \"msLevel 2\" --filter \"zeroSample removeExtra\" --outdir \"#{local_base}/out/#{bacterium}\"";
		puts `"#{msconvert}" #{raw} --mgf --filter "msLevel 2" --filter "zeroSample removeExtra" --outdir "#{local_base}/out/#{bacterium}"`;


		mgf = raw.gsub(/RAW$/,"mgf")

		File.delete(raw)

		ftp.chdir("#{remote_base}/sequence/#{bacterium}")
		ftp.nlst.each do |sequence|

			FileUtils.mkdir_p "#{local_base}/sequence/#{bacterium}" unless File.directory? "#{local_base}/sequence/#{bacterium}"
			Dir::chdir "#{local_base}/sequence/#{bacterium}"

			ftp.chdir("#{remote_base}/sequence/#{bacterium}")

			if File.exists? sequence
				puts "--#{sequence} already exists"
			else
				puts "--Getting #{sequence}"
				ftp.get(sequence)
			end


			outdir = "#{local_base}/out/#{bacterium}/#{sequence.gsub(/fasta$/,'')}"
			FileUtils.mkdir_p "#{outdir}" unless File.directory? "#{outdir}"
			Dir::chdir "#{outdir}"

			puts "--Creating config"
			config = mgf.gsub(/mgf$/,"config")
			File.open(config,"w") do |f|

				f.puts "Spectra=#{local_base}/out/#{bacterium}/#{mgf}"
				f.puts "Fasta=#{local_base}/sequence/#{bacterium}/#{sequence}"

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




			File.delete(sequence)

		end	#	ftp.nlst.each do |sequence|

	end	#	ftp.nlst('*').each |raw|

end	#	ftp.nlst('*').each do |bacterium|


ftp.close
puts "Closed"
