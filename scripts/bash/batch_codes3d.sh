#!/bin/bash

snp_file="../../data/all_gwas_snps.txt"
snpdir="../../data/batch_snps"
#mkdir $snpdir
#cat $snp_file | (cd $snpdir; split -l 50 -d)
results_dir="../../results/batch_results"
#mkdir -p $results_dir
done_dir="../../data/batch_snps/done"
#mkdir $done_dir
for f in $snpdir/*
do
    if [[ -f "${f}" ]]
    then
	echo Processing $f 
	python ../python/codes3d/codes3d/codes3d.py -i $f\
	       -o $results_dir/${f:${#snpdir}+1} -r mboi
	mv $snpdir/${f:${#snpdir}+1} $done_dir/
    fi
done
