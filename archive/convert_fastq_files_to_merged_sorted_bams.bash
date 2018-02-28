#!/usr/bin/env bash


for fastq1 in $( ls *R1_001.fastq.gz ) ; do

	basename=${fastq1%%_R1_001.fastq.gz}
	fastq2=${basename}_R2_001.fastq.gz

	java -jar ~/picard.jar FastqToSam \
		F1=$fastq1 \
		F2=$fastq2 \
		O=$basename.bam \
		SM=$basename

done


#	ASSUME_SORTED=true \	#	What if they aren't???
#	Guessing that if assume sorted, doesn't try to sort them.
#	Does drop seq expect them to be sorted?

java -jar ~/picard.jar MergeSamFiles \
	ASSUME_SORTED=false \
	SORT_ORDER=queryname \
	OUTPUT=$PWD.bam \
	$( ls *bam | awk '{printf " INPUT="$0}' )

