#!/usr/bin/env bash


calling_dir=$PWD


while [ $# -ne 0 ] ; do
	cd $calling_dir
	echo "Processing :${1}:"

	mkdir -p ${1}
	cd ${1}

	cmd="seurat_group.R $1"
	echo $cmd
	$cmd

	shift
done
