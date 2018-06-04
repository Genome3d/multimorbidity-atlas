#!/usr/bin/python

import csv
import os
import random
import argparse
"""Computes matrices of shared eQTLs and eGene among traits, as well
    as null gene matrix.

    Output:
        egene_ratios.txt: A matrix of shared eGenes among phenotypes
        eqtl_ratios.txt: A matrix of shared eQTLs among phenotypes
        spatial_phenotypes.txt: eQTL associations of phenotypes that do not
            share eQTLs but share target eGenes.
"""

def parse_data(data_fp):
    """Parse significat eQTL association files for traits.

    Args:
        data_fp: Filepath of directory containing trait directories, which in
            turn contain significant_eqtls.txt files from CoDeS3D.

    Output:
        traits: A dictionary of traits, their gene and eQTL lists
            {'acne':
                {'genes': ['ABO', 'FOX2'],
                 'gene_number': 2,
                 'snps': ['rs1234'],
                 'snp_number': 1},
            'obesity:
                {...}
            }

        all_genes: A list of all eGenes from all traits in the input directory.
    """
    traits = {}
    eqtls = {}
    all_genes = []
    _,traits_dir,_ = next(os.walk(data_fp), (None,[],None))
    traits_to_remove = []
    print('Parsing eQTL data for..')
    for trait in traits_dir:
        #trait = t[:len(t)-4]
        print('\t'+trait)
        if trait not in traits.keys():
            traits[trait] = {'genes':[], 'gene_number':0, 'snps':[], 'snp_number':0}
        fp = open(os.path.join(data_fp,trait, 'significant_eqtls.txt'), 'r')
        reader = csv.reader(fp, delimiter = '\t')
        next(reader, None)
        for row in reader:
            gene = row[3]
            snp = row[0]
            if gene not in traits[trait]['genes']:
                traits[trait]['genes'].append(gene)
                traits[trait]['gene_number'] += 1
            if snp not in traits[trait]['snps']:
                traits[trait]['snps'].append(snp)
                traits[trait]['snp_number'] += 1
            if gene not in all_genes:
                all_genes.append(gene)
        if not traits[trait]['genes'] or not traits[trait]['snps']:
            traits_to_remove.append(trait)
    for trait in traits_to_remove:
        del traits[trait]  #Omit traits with no siginicant eQTL associations
    return traits, all_genes

def calc_ratios(traits):
    """Compute matrix of shared eQTLs and eGene among traits.

    Args:
        traits: A dictionary from parse_data.

    Output:
        egene_ratios.txt: A matrix of shared eGenes among phenotypes
        eqtl_ratios.txt: A matrix of shared eQTLs among phenotypes
        spatial_phenotypes.txt: eQTL associations of phenotypes that do not
            share eQTLs but share target eGenes.
    """
    print('Calculating ratios..')
    outfile = open(os.path.join(args.output,'egene_matrix.txt'), 'w')
    gene_writer = csv.writer(outfile, delimiter = '\t')
    gene_writer.writerow(['Trait1', 'Trait2', 'Genes1', 'Genes2', 'Ratio'])
    snpfile = open(os.path.join(args.output, 'eqtl_matrix.txt'), 'w')
    snp_writer = csv.writer(snpfile, delimiter = '\t')
    snp_writer.writerow(['Trait1', 'Trait2', 'eQTLs1', 'eQTLs2', 'Ratio'])
    spatial_file = open(os.path.join(args.output, 'spatial_traits.txt'), 'w')
    spatial_writer = csv.writer(spatial_file, delimiter = '\t')
    spatial_writer.writerow(['Trait1', 'Trait2', 'Common_eGenes', 'eGene_Ratio',
                             'eGene_List'])
    trait_file = open(os.path.join(args.output, 'matrix_traits.txt'), 'w')
    trait_file.write('Traits\n')

    for trait1 in sorted(traits.keys()):
        trait_file.write(trait1 + '\n')
        for trait2 in sorted(traits.keys()):
            common_genes = list(set(traits[trait1]['genes']) &
                                set(traits[trait2]['genes']))
            gene_ratio = float(len(common_genes)) / traits[trait1]['gene_number']
            gene_writer.writerow((trait1, trait2, traits[trait1]['gene_number'],
                                  traits[trait2]['gene_number'],
                                  len(common_genes), gene_ratio))
            common_snps = list(set(traits[trait1]['snps']) &
                               set(traits[trait2]['snps']))
            snp_ratio = float(len(common_snps)) / traits[trait1]['snp_number']
            snp_writer.writerow((trait1, trait2, traits[trait1]['snp_number'],
                             traits[trait2]['snp_number'], len(common_snps),
                                 snp_ratio))
            if snp_ratio == 0 and gene_ratio > 0:
                spatial_writer.writerow((trait1, trait2, len(common_genes),
                                         gene_ratio, common_genes))
    outfile.close()
    snpfile.close()
    spatial_file.close()
    trait_file.close()

def simulate_nulls(traits, gene_pool):
    """Generate 100 null sets of eGene ratio matrix.

    Args:
        traits: Dicitionary of traits and their lists of eGenes and eQTLs.
        gene_pool: A list of all eGenes associated with traits being analysed.

    Output:
        control_gene_ratios.txt: A matrix of the mean eGene ratios of the 100
            null matrices generated.
    """
    replicates = {}
    print('Simulating null dataset..')
    null_sets = 100
    for i in range(null_sets):
        print('\t'+ str(i))
        rand_traits = {}
        for trait in traits:
            if trait not in rand_traits:
                rand_traits[trait] = {'genes':[], 'gene_number':0}
            j = 0
            while j  <= traits[trait]['gene_number']:
                gene = random.choice(gene_pool)
                if gene not in rand_traits[trait]['genes']:
                    rand_traits[trait]['genes'].append(gene)
                    rand_traits[trait]['gene_number'] += 1
                    j += 1
        for trait1 in rand_traits:
            if trait1 not in replicates:
                replicates[trait1] = {}
            for trait2 in rand_traits:
                common = list(set(rand_traits[trait1]['genes']) & set(rand_traits[trait2]['genes']))
                ratio = float(len(common)) / rand_traits[trait1]['gene_number']
                if trait2 not in replicates[trait1]:
                    replicates[trait1][trait2] = 0
                replicates[trait1][trait2] +=ratio
    outfile = open(os.path.join(args.output, 'control_gene_ratios.txt'), 'w')
    writer = csv.writer(outfile, delimiter = '\t')
    writer.writerow(['Trait1', 'Trait2', 'Ratio'])
    for trait1 in replicates:
        for trait2 in replicates[trait1]:
            writer.writerow((trait1, trait2, replicates[trait1][trait2]/null_sets))
    outfile.close()
    print('Simulation completed.')


if __name__=='__main__':
    parser = argparse.ArgumentParser(description = '')
    parser.add_argument('-i', '--input', default='../../results/trait_eqtls',
                        help='Directory of trait directories containing' +\
                        'significant spatial eQTL-eGene associations.')
    parser.add_argument('-o', '--output', default='../../results/',
                        help='Filepath for output')
    args = parser.parse_args()
    traits, all_genes = parse_data(args.input)
    calc_ratios(traits)
    #simulate_nulls(traits, all_genes)
