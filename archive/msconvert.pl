#!/usr/bin/env perl

use strict;
use warnings;
use Net::FTP;
use File::Path qw(make_path remove_tree);

my $ftp = Net::FTP->new("massive.ucsd.edu", Debug => 0)
  or die "Cannot connect to host: $@";

#	Apparently NEED to login when scripting
$ftp->login("anonymous",'-anonymous@')
	or die "Cannot login ", $ftp->message;

my $local_root  = "/Users/jakewendt/";
my $remote_root = "/";
my $base = "MSV000079053";
my $local_base  = $local_root.$base;
my $remote_base = $remote_root.$base;

$ftp->cwd($remote_base."/raw")
	or die "Cannot change working directory ", $ftp->message;
my @bacteria = $ftp->ls;

if( ! -d $local_base."/raw" ) {
	make_path( $local_base."/raw" ) or die "Failed make_path $!";
}

if( ! -d $local_base."/sequence" ){
	make_path( $local_base."/sequence" ) or die "Failed make_path $!";
}

if( ! -d $local_base."/out" ){
	make_path( $local_base."/out" ) or die "Failed make_path $!";
}


#	http://fields.scripps.edu/rawconv/download/MSFileReader%202.2.62.zip

my $msconvert = 'C:\Program Files\ProteoWizard\ProteoWizard 3.0.18187.b51377ef8\msconvert.exe';

foreach (@bacteria) {
	my $bacterium = $_;
	print "$bacterium\n";

	$ftp->cwd($remote_base."/raw/".$bacterium)
		or die "Cannot change working directory ", $ftp->message;
	my @raws = $ftp->ls;

	foreach(@raws) {
		my $raw = $_;
		print "$raw\n";
		chdir $local_base."/raw/";
		mkdir $bacterium if( ! -d $bacterium );
		chdir $bacterium;

		$ftp->cwd($remote_base."/raw/".$bacterium)
			or die "Cannot change working directory ", $ftp->message;
		if( ! -e $raw ){
			$ftp->get($raw)
				or die "get raw file failed ", $ftp->message;
		}

		#	lots of quotes are needed
		print `"$msconvert" $raw --mgf --filter "msLevel 2" --filter "zeroSample removeExtra"`;


		$ftp->cwd($remote_base."/sequence/".$bacterium)
			or die "Cannot change working directory ", $ftp->message;

		my @sequences = $ftp->ls;
#		die "Not just 1 sequence found." if $#sequences != 0;
#		Some have more than 1 sequence.
#		die "No sequence found." if $#sequences < 0;

		foreach(@sequences) {
			my $sequence = $_;
			print "$sequence\n";

			chdir $local_base."/sequence/";
			mkdir $bacterium if( ! -d $bacterium );
			chdir $bacterium;

			$ftp->cwd($remote_base."/sequence/".$bacterium)
				or die "Cannot change working directory ", $ftp->message;
			if( ! -e $sequence ){
				$ftp->get($sequence)
					or die "get sequence failed ", $ftp->message;
			}


			#	MODa

#print `ls -ali`;



		}	#	foreach(@sequences)

	}	#	foreach(@raws) {
}

$ftp->quit;

