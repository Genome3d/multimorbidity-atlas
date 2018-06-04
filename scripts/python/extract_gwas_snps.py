#! /usr/bin/env python

import sys
import csv
import os
import argparse


def get_catalogue_snps(gwas_fp):
    """
    Extract SNPs from all associations downloaded from GWAS Catalogue.
    """
    trait_snps = {}
    all_snps = []
    snps_bed = []
    error_log = []
    if os.path.isfile(gwas_fp):
        with open(gwas_fp, 'r') as gwas_file:
            associations = csv.reader(gwas_file, delimiter = '\t') 
            next(associations, None)
            for line in associations:
                snps = []
                chromosomes = []
                snp_locus_column = line[12].split(';')
                snps_column = line[21]
                trait = line[34].replace(' ', '_')
                trait = trait.replace('/', '-')
                trait = trait.replace(',', '_')
                trait = trait.replace(':', '_')
                
                chr_column = line[11].split(';')
                if chr_column is not None:
                    for chrom in chr_column:
                        snp_chr = 'chr' + chrom
                        chromosomes.append(snp_chr)
                if line[21].startswith('rs'):
                    snps_column = snps_column.strip().split(';')
                    snp_x = []
                    chrom_x = []
                    locus_x = []
                    for snp in snps_column:
                        if len(snp.split()) > 1: 
                           snp_x.append(snp.replace('x', ''))
                    if snp_x:
                        snps_column = snp_x
                    for chro in chromosomes:
                        if len(chro.split()) > 1: 
                            chro = chro.replace('x', '')
                            chs = []
                            for ch in chro.split():
                                if not ch.startswith('chr'):
                                    ch = 'chr' + ch
                                    chs.append(ch)
                                else:
                                    chs.append(ch)
                            chro = chs
                            chrom_x = chro
                    if chrom_x:
                        chromosomes = chrom_x
                    for loc in snp_locus_column:
                        if len(loc.split()) > 1: 
                           loc = loc.replace('x', ' ').strip()
                           locus_x = loc.split()
                    if locus_x:
                        snp_locus_column = locus_x

                    for i in range(0, len(snps_column)):
                        snp = snps_column[i].strip()
                        snp_chr = ''
                        snp_locus = ''
                        if len(chromosomes) == 1 and len(snps_column) == 1:
                            snp_chr = chromosomes[0]
                            snp_locus = snp_locus_column[0]
                        elif len(chromosomes) == len(snps_column):
                            snp_chr = chromosomes[i]
                            snp_locus = snp_locus_column[i]
                        elif len(chromosomes) == len(snp.split()) \
                             and len(snp.split()) > 1:
                            split_snps = snp.split()
                            for r in range(0, len(split_snps)):
                                if r == 0:
                                    snp = split_snps[r]
                                    snp_chr = chromosomes[i]
                                    snp_locus = snp_locus_column[i]
                                else:
                                    snp = split_snps[r]
                                    snp_chr = chromosomes[r]
                                    snp_locus = snp_locus_column[r]
                                if (not snp in snps) and snp_locus != '':
                                    snps.append(snp)
                                    to_trait_snps = []
                                    if trait not in trait_snps.keys():
                                        to_trait_snps.append(snp)
                                        trait_snps[trait] = to_trait_snps
                                    else:
                                        to_trait_snps = trait_snps[trait]
                                        to_trait_snps.append(snp)
                                        trait_snps[trait] = to_trait_snps
                                if snp not in all_snps and snp_locus != '':
                                    all_snps.append(snp.strip())
                                    to_bed = snp_chr, int(snp_locus)-1, snp_locus, snp
                                    snps_bed.append(to_bed)
                                continue
                        elif len(snp.split()) == 1 and len(chromosomes) == 2:
                            if chromosomes[0] == chromosomes[1]:
                                snp_chr = chromosomes[1]
                                loci_0 = int(snp_locus_column[0])
                                loci_1 = int(snp_locus_column[1])
                                if (loci_1 - loci_0) == 1:
                                    snp_locus = loci_1
                                else:
                                    error_log.append(trait, snp, snp_chr, snp_locus)
                            else:
                                error_log.append(trait, snp, snp_chr, snp_locus)
                        else:
                            # TODO: fix uncorresponding SNPs and chromosomes in a trait.
                            rectify = trait, snp, snp_chr, snp_locus
                            error_log.append(rectify)
                        if (not snp in snps) and snp_locus != '':
                            snps.append(snp)
                            to_trait_snps = []
                            if trait not in trait_snps.keys():
                                to_trait_snps.append(snp)
                                trait_snps[trait] = to_trait_snps
                            else:
                                to_trait_snps = trait_snps[trait]
                                to_trait_snps.append(snp)
                                trait_snps[trait] = to_trait_snps
                        if snp not in all_snps and snp_locus != '':
                            all_snps.append(snp.strip())
                            to_bed = snp_chr, int(snp_locus)-1, snp_locus, snp
                            snps_bed.append(to_bed)
    else:
        print(gwas_fp + ' not a  file.')

    trait_snps_dir = os.path.join(output_fp, 'trait_snps')
    if not os.path.isdir(trait_snps_dir):
        os.mkdir(trait_snps_dir)
    for trait in trait_snps:
        trait_file = open (os.path.join(trait_snps_dir, trait + '.txt'), 'w')
        for snp in trait_snps[trait]:
            trait_file.write(snp + '\n')
            if not snp in all_snps:
                all_snps.append(snp)
        trait_file.close()
    allsnps = open(os.path.join(output_fp, 'all_gwas_snps.txt'), 'w')
    for snp in all_snps:
        allsnps.write(snp + '\n')
    allsnps.close()
    with open(output_fp + '/all_gwas_snps.bed', 'w') as bed_file:
        writer = csv.writer(bed_file, delimiter = '\t')
        writer.writerow(['chrom', 'chromStart', 'chromEnd', 'name'])
        writer.writerows(snps_bed)
    bed_file.close()
    with open(output_fp + '/error_log.txt', 'w') as error_file:
        writer = csv.writer(error_file, delimiter = '\t')
        writer.writerow(['trait', 'snp'])
        writer.writerows(error_log)
    error_file.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "")
    parser.add_argument("-i", "--input",
                        default='../../data/downloads/' +\
                        'gwas_catalog_v1.0.1-associations_e85_r2016-08-25.tsv',
                        help = "Filepath of directory containing downloaded "+\
                        "GWAS association from GWAS Catalogue.")
    parser.add_argument("-o", "--output",
                        default='../../data/',
                        help = "Filepath of output directory.")
    args = parser.parse_args()
    gwas_assoc_fp = args.input
    output_fp = args.output
    get_catalogue_snps(gwas_assoc_fp)
