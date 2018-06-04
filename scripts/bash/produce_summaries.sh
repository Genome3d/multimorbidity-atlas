#!/bin/bash
<<'COMMENT'
Computes adjusted p values using the Benjamini-Hochberg FDR and generates summary
files for each trait in ../../results/trait_eqtls/

Args: eqtls.txt file for each trait

Output: summary.txt file containing spatial eQTL-gene interactions in tissues.
COMMENT
repdir=../../results/trait_eqtls
echo 'Producing summary for...'
for i in $(ls $repdir);
do
    echo -e '\t'$i
    python ../python/codes3d/codes3d/produce_summary.py -e $repdir/$i/eqtls.txt -o $repdir/$i \
	   -c ../python/codes3d/docs/codes3d.conf
done
