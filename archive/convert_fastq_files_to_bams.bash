#!/usr/bin/env bash

basedir=$PWD
for dir in $( ls -d 1763* ) ; do
	cd $basedir

	echo $dir
	cd $dir
#	ls *fastq*
#2B_S1_L004_R1_001.fastq.gz  2B_S1_L004_R2_001.fastq.gz

	fastq1=$( ls *R1*fastq.gz )
	fastq2=$( ls *R2*fastq.gz )
	basename=${fastq1%%_R1_001.fastq.gz}
	echo $basename

	java -jar $basedir/picard.jar FastqToSam \
		F1=$fastq1 \
		F2=$fastq2 \
		O=$basename.bam \
		SM=$basename

done
