#!/usr/bin/python

import os
import csv
import argparse
import sys
import sqlite3
"""
Get multimorbid traits (and their eGenes) associated with the query trait
"""

def parse_traits(traits_fp):
    """Returns list of traits in spatial eGene matrix.

    Args:
        traits_fp: ../../results/matrix_file.txt
    """
    traits = []
    with open(traits_fp, 'r') as traits_file:
        traits = [t for t in (line.strip() for line in traits_file) if t]
    return traits
    
def read_matrix(matrix_file):
    """Parses shared eGene ratios.

    Args:
        matrix_file: ../../results/egene_matrix.txt

    Returns:
        matrix: A list of rows from the input file.
    """
    matrix = []
    mfile = open(matrix_file, 'r')
    mreader = csv.reader(mfile, delimiter = '\t')
    next(mreader, None)
    for row in mreader:
        matrix.append(row)
        if not row[0] in matrix_traits:
                matrix_traits.append(row[0])
        if not row[1] in matrix_traits:
            matrix_traits.append(row[1])
    return matrix

def parse_input(to_query, matrix_traits):
    """Parses and validates user input.

    Args:
        to_query: An text file containing a list of traits, SNPs or genes to
            query. Or a space-delimited list of same from the terminal with -i

        matrix_traits: A list of traits in the eGene matrix.

    Returns:
        query_traits: A validated list of traits that are associated with user
            input.
    """
    conn = sqlite3.connect(args.eqtlsDB)
    conn.text_factory = str
    cur = conn.cursor()
    query_input = []
    query_traits = []
    if os.path.isfile(to_query[0]):  # Read traits from file
        with open(to_query[0], 'r') as input_file:
            query_input = [t for t in (line.strip() for line in input_file) if t]
    else:  # Read space-delimited traits
        query_input = to_query
    upper_matrix_traits = [t.upper() for t in matrix_traits]
    for item in query_input:
        if item.startswith('rs'):
            cur.execute("SELECT DISTINCT trait FROM interactions WHERE snp = ?;",
                        (item,))
            from_cur = [trait[0] for trait in cur.fetchall()]
            query_traits += from_cur
        else:
            cur.execute("SELECT DISTINCT trait FROM interactions WHERE gene = ?;",
                        (item,))
            from_cur = [trait[0] for trait in cur.fetchall()]
            if from_cur:
                query_traits += from_cur
            else:
                try:
                    idx = upper_matrix_traits.index(item.upper())
                    query_traits.append(matrix_traits[idx])
                except:
                    query_traits.append(item)
    if not set(query_traits) <= set(matrix_traits):
        #Handle user input error
        sys.exit('\nOops! Program is terminated because the following traits are '+\
                 'not found in the gene matrix.\n'+
                 'Please use only traits in "../../results/matrix_traits.txt"\n\t'+
                 '\n\t'.join(trait for trait in query_traits if not trait in
                             matrix_traits)+'\n')
    return query_traits
    
def find_traits(query_traits, matrix):
    """Identifies traits that share spatial eGenes and gives a score to each pair.

    Args:
        query_traits: Validated list of traits from user input.
        matrix: A list of trait-pairs that share eGenes.

    Returns:
        results:
        snp_data_c:
    """
    comorbid = []
    traits = {}
    results = []
    snp_data = []
    snp_data_c = []
    gene_bank = {}
    for trait in query_traits:
        for row in matrix:
            if (row[0] in query_traits) and (row[1] in query_traits):
                comorbid.append(row)
    for row in comorbid:
        xtrait = row[0]
        ytrait = row[1]
        if xtrait not in traits:
            traits[xtrait] = get_genes(xtrait)
        if ytrait not in traits:
            traits[ytrait] = get_genes(ytrait)
        common_genes = list(set(traits[xtrait].keys()).intersection(traits[ytrait].keys()))
        for gene in common_genes:
            if gene not in gene_bank.keys():
                gene_bank[gene]={'count':0,'snps':{}}
    for egene in gene_bank:
        for trait in traits:
            if egene in traits[trait].keys():
                gene_bank[egene]['count'] += 1
                snps = traits[trait][egene]['snps']
                for snp in snps:
                    if snp not in gene_bank[egene]['snps'].keys():
                        gene_bank[egene]['snps'][snp]=[]
                    gene_bank[egene]['snps'][snp].append(trait)
    for egene in gene_bank:
        gene_score = gene_bank[egene]['count']/float(len(query_traits))
        snp_list = ''
        i=0
        for snp in gene_bank[egene]['snps']:
            i +=1
            traits = gene_bank[egene]['snps'][snp]
            snp_data.append(get_snps(egene,gene_bank[egene]['snps']))
            trait_list = traits[0]
            snp_list = snp_list + snp + ': '
            if len(traits) > 1:
                for i in range(1,len(traits)):
                    trait_list = trait_list + ', ' + traits[i]
            if i < len(gene_bank[egene]['snps']):
                trait_list = trait_list + '; '
            snp_list = snp_list + trait_list
        results.append((egene, gene_bank[egene]['count'], "{0:.2g}".format(gene_score), snp_list))
    for row in snp_data:
        if row not in snp_data_c:
            snp_data_c.append(row)
    return(results, snp_data_c)


def get_snps(egene, egene_dict):
    snp_data = []
    for eqtl in egene_dict:
        traits = egene_dict[eqtl]
        trait_list = traits[0]
        if len(trait_list) > 1:
            for i in range(1, len(traits)):
                trait_list += ', ' + traits[i]
        for trait in traits:
            trait_dir = os.path.join(args.traits_dir, trait)
            eqtl_file = open(os.path.join(trait_dir, 'significant_eqtls.txt'), 'r')
            ereader = csv.reader(eqtl_file, delimiter = '\t')
            next(ereader, None)
            for row in ereader:
                if row[0] == eqtl and row[3] == egene:
                    gene = row[3]
                    #gene_chr = row[4]
                    #gene_start = row[5]
                    #gene_end = row[6]
                    snp = row[0]
                    snp_chr = row[1]
                    snp_locus = row[2]
                    distance = row[14]
                    to_ = [snp, snp_chr, snp_locus, gene, distance, trait_list]
                    if to_ not in snp_data:
                        snp_data.append(to_)
                    break
            break
    return(snp_data)

def get_genes(trait):
    trait_dict = {}
    trait_dir = os.path.join(args.traits_dir, trait)
    eqtl_file = open(os.path.join(trait_dir, 'significant_eqtls.txt'), 'r')
    ereader = csv.reader(eqtl_file, delimiter = '\t')
    next(ereader, None)
    for row in ereader:
        gene = row[3]
        snp = row[0]
        tissue = row[7]
        if gene not in trait_dict:
            trait_dict[gene] = {'tissues':[], 'snps':[]}
        if snp not in trait_dict[gene]['snps']:
            trait_dict[gene]['snps'].append(snp)
        if tissue not in trait_dict[gene]['tissues']:
            trait_dict[gene]['tissues'].append(tissue)
    return(trait_dict)

def write_to_file(results, snp_data, outputdir):
    if not os.path.isdir(outputdir):
        os.mkdir(outputdir)
    rfile = open(os.path.join(outputdir, 'cluster_analysis.txt'), 'w')
    rwriter = csv.writer(rfile, delimiter = '\t')
    rwriter.writerow(['eGene', 'Traits', 'eGene_Score', 'eQTLs'])
    rwriter.writerows(results)
    rfile.close()

    sfile = open(os.path.join(outputdir, 'eQTL_interactions.txt'), 'w')
    swriter = csv.writer(sfile, delimiter = '\t')
    swriter.writerow(['eQTL', 'eQTL_Chr', 'eQTL_Locus', 'eGene', 'Distance', 'Traits'])
    for row in snp_data:
        swriter.writerows(row)
    sfile.close()


if __name__=='__main__':
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('-i', '--input', nargs='+', required = True, \
                            help = 'Traits to be queried. Please refer to '+\
                        '\"../../results/matrix_traits.txt\" for list of '+\
                        'traits. SNP rsIds and gene names are also accepted.')
    parser.add_argument('-m', '--matrix_file',
                        default='../../results/egene_matrix.txt',
                        help='\"egene_matrix.txt\", the output file from '+\
                        '\"./control_analysis.py\"')
    parser.add_argument('-t', '--matrix_traits',
                        default='../../results/matrix_traits.txt',
                        help='\"matrix_traits.txt\", containing list of '+\
                        'traits in eGene matrix.')
    parser.add_argument('-s', '--traits_dir',
                        default='../../../results/trait_eqtls/',
                        help='Filepath to trait directories with '+
                        '\"significant_eqtls.txt\" files.')
    parser.add_argument('-e', '--eqtlsDB',
                        default='../../results/spatial_eqtls.db',
                        help='Spatial eQTL database \"spatial_eqtls.db\".')
    parser.add_argument('-o', '--output', required = True, \
                            help = 'Filepath of output directory')
    args = parser.parse_args()
    matrix_traits = parse_traits(args.matrix_traits)
    query_traits = parse_input(args.input, matrix_traits)
    matrix = read_matrix(args.matrix_file)
    results, snp_data = find_traits(query_traits, matrix)
    write_to_file(results,snp_data,args.output)
