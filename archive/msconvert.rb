#!/usr/bin/env ruby

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



#
##	http://fields.scripps.edu/rawconv/download/MSFileReader%202.2.62.zip
#
#my $msconvert = 'C:\Program Files\ProteoWizard\ProteoWizard 3.0.18187.b51377ef8\msconvert.exe';
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

#		#	lots of quotes are needed
#		print `"$msconvert" $raw --mgf --filter "msLevel 2" --filter "zeroSample removeExtra"`;






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





			puts "--Creating config"


			puts "--Running MODa on #{raw} and #{sequence}"




		end	#	ftp.nlst.each do |sequence|

	end	#	ftp.nlst('*').each |raw|

end	#	ftp.nlst('*').each do |bacterium|


ftp.close
puts "Closed"
