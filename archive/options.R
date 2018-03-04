#!/usr/bin/env Rscript

#	install.packages('optparse')
#	http://tuxette.nathalievilla.org/?p=1696
library("optparse")
 
#	default action is "store"
#	default dest variable are the options without the leading dashes
option_list = list(
	make_option(c("-m", "--manhattan"), type="character", default=NULL, #	dest= m and manhattan
		help="manhattan.plot file name", metavar="character"),
	make_option(c("-q", "--qq"), type="character", default=NULL, action="store", #dest="xyz", = q and qq
		help="qq file name", metavar="character"),

  make_option(c("--quiet"), action="store_true", default=FALSE,
              help="Make the program not be verbose."),

  make_option("--test", action="store_true", default=FALSE,
              help="testing."),

	make_option(c("-o", "--outpath"), type="character", default="./", 
		help="output file name [default= %default]", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

print(opt$manhattan)
print(opt$test)
print(opt$q)
print(opt$qq)
print(opt$xyz)
print(opt$quiet)

