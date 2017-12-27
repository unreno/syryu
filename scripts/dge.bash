#!/usr/bin/env bash



#	find . -name error_detected.bam -execdir ~/syryu/drop_seq_alignment/dge.bash \; > dge.log 2>&1


echo
echo
echo
pwd
echo



~/singlecell/Drop-seq_tools-1.13/BAMTagHistogram \
	INPUT=error_detected.bam \
	OUTPUT=out_cell_readcounts.txt.gz \
	TAG=XC

#
#	Only keep those with more than one. This is just a test.
#
#zcat out_cell_readcounts.txt.gz | tail -n +2 | awk '( $1 > 1 ){print $2}' | gzip > cell_bc_file.txt.gz
zcat out_cell_readcounts.txt.gz | tail -n +2 | awk '{print $2}' | gzip > cell_bc_file.txt.gz

~/singlecell/Drop-seq_tools-1.13/DigitalExpression \
	INPUT=error_detected.bam \
	OUTPUT=error_detected.dge.txt.gz \
	CELL_BC_FILE=cell_bc_file.txt.gz \
	SUMMARY=out_gene_exon_tagged.dge.summary.txt \
	NUM_CORE_BARCODES=100


#	NUM_CORE_BARCODES=100 ... Doesn't seem to make any difference. (Maybe default value?)

./seurat.R

echo
