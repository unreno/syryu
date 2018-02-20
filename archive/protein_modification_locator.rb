#!/usr/bin/env ruby

require 'csv'

out=File.open("ProteinModification.txt",'w')
out.puts "Leading razor protein\tAA\tAA location\tSequences with Modified AA\tSequences with Unmodified AA" 

CSV.open("MatchedModification.txt",'rb',
	{ headers: true, col_sep: "\t"  }).each do |line|

	puts line.to_hash





end






out.close
