#!/usr/bin/python
import os
import csv
import sqlite3
import sys
import argparse

"""Creates eQTL association table in ../../results/spatial_eqtls.db
   Writes gene associations to ../../results/genes_details.txt
"""
def parse_files(traits_dir, gwas_genes):
    """
    
    Args:
        traits_dir: ../../results/trait_eqtls/
        gwas_genes: A dictionary of genes mapped to SNPs in the GWAS Catalog.

    Output:
        interactions: A table in ../../results/spatial_eqtls.db with ff schema:
            1. trait  # Name of (composite) GWAS trait for which SNP is mapped.
            2. gene   # Spatially mapped eGene.
            3. snp    # Spatial eQTL.
            4. cis    # False if eQTL-eGene are >1MB or inter-chromosomal. 
            5. mapped # 0 if snp is not mapped to gene in the GWAS Catalog.
    """
    all_eqtls = []
    gene_checked = []
    eqtls = {}
    egenes = {}
    eqtl_list = []
    egene_list = {}
    conn = sqlite3.connect(args.eqtlsDB)
    conn.text_factory == str
    cur = conn.cursor()
    cur.execute("DROP TABLE IF EXISTS interactions;")
    cur.execute("""
    CREATE TABLE interactions (trait TEXT, gene TEXT, snp TEXT, 
    cis TEXT, mapped TEXT)
    """)
    print('Extracting eQTL interactions in...')
    for trait in os.listdir(traits_dir):
        eqtl_fp = os.path.join(traits_dir, trait, 'significant_eqtls.txt')
        if os.path.isfile(eqtl_fp):
            eqtlfile = open(eqtl_fp, 'r')
            freader = csv.reader(eqtlfile, delimiter = '\t')
            next(freader, None)
            done_list = []
            print('\t' + trait)
            for row in freader:
                mapped = False
                if row[3] in gwas_genes.keys():
                    for snp in gwas_genes[row[3]]['snps']:
                        if row[0] == snp[0]:
                            mapped = True
                            break
                snpchr = row[1]
                gene = row[3]
                snp = row[0]
                tissue = row[7]
                interaction = row[13]
                if gene not in egene_list.keys():
                    egene_list[gene] = {'snps':[], 'traits':[], 'tissues':[], 'interactions':[]}
                if snp not in egene_list[gene]['snps']:
                    egene_list[gene]['snps'].append(snp)
                if trait not in egene_list[gene]['traits']:
                    egene_list[gene]['traits'].append(trait)
                if tissue  not in egene_list[gene]['tissues']:
                    egene_list[gene]['tissues'].append(tissue)
                if interaction not in egene_list[gene]['interactions']:
                    egene_list[gene]['interactions'].append(interaction)
                if mapped:
                    egene_list[gene]['mapped'] = 'True'
                else:
                    egene_list[gene]['mapped'] = 'False'
                if not snp+gene in done_list:
                    cur.execute("INSERT INTO interactions VALUES (?,?,?,?,?);",
                                (trait, gene, snp, interaction, mapped))
                    done_list.append(snp+gene)
    conn.commit()
    conn.close()
    rfile = open(os.path.join(args.output, 'genes_details.txt'), 'w')
    rwriter = csv.writer(rfile, delimiter= '\t')
    rwriter.writerow(('Gene', 'SNPs', 'Traits', 'Tissues', 'Interactions',
                      'Mapped'))
    for gene in egene_list:
        interactions = ''
        inter = egene_list[gene]['interactions']
        if len(inter) > 1:
            interactions = 'both'
        elif inter[0] == 'True':
            interactions = 'cis'
        else:
            interactions = 'trans'
        snplist = egene_list[gene]['snps']
        snps = snplist[0]
        if len(snplist) > 1:
            for i in range(1, len(snplist)):
                snps = snps + ', ' + snplist[i]
        traitlist = egene_list[gene]['traits']
        traits = traitlist[0]
        if len(traitlist) > 1:
            for i in range(1, len(traitlist)):
                traits = traits + ', ' + traitlist[i]
        tissuelist = egene_list[gene]['tissues']
        tissues = tissuelist[0]
        if len(tissuelist) > 1:
            for i in range(1, len(tissuelist)):
                tissues = tissues + ', ' + tissuelist[i]
        rwriter.writerow((gene, snps, traits, tissues,interactions,
                          egene_list[gene]['mapped']))
    rfile.close()


def parse_gwas(gwas_fp):
    """Extract mapped genes in GWAS Catalog associations

    Args:
        gwas_fp: Filepath to GWAS associations file.

    Returns:
        gwas_genes: A dictionary with ff structure:
            {
             'ABO': {
                'snps': ['rs1234', 'rs2345'],
                'traits': ['acne', 'type_II_diabetes_mellitus']
                }
             'FTO': {...}
             }
    """
    gwas_genes = {}
    testers = []
    gwas_file = open(gwas_fp, 'r')
    reader = csv.reader(gwas_file, delimiter = '\t')
    next(reader, None)
    print('Parsing GWAS Catalog associations...')
    for row in reader:
        snp = row[21].strip()
        snp = snp.split(';')
        genes = row[14]
        genes = row[14].replace(';',',')
        genes = genes.split(',')
        trait = row[34]
        container = []
        for gene in genes:
            gene = gene.strip()
            if gene != 'No mapped genes' and gene != '':
                gene = gene.strip()
                tester = gene+str(snp)+trait
                if not tester in testers:
                    testers.append(tester)
                    if gene not in gwas_genes.keys():
                        gwas_genes[gene] = {'snps':[],'trait':[]}
                    if snp not in gwas_genes[gene]['snps']:
                        gwas_genes[gene]['snps'].append(snp)
                    if trait not in gwas_genes[gene]['trait']:
                        gwas_genes[gene]['trait'].append(trait)
    #print(len(gwas_genes))
    return(gwas_genes)

def parse_hic(hic, eqtls):
    
    hic_file=open(hic, 'r')
    reader = csv.reader(hic_file, delimiter='\t')
    next(reader, None)
    hic_done = []
    to_file = []
    for row in reader:
        eqtl_egene = row[0]+'_'+row[1]
        eqtl_done = []
        cells = row[3].split(',')
        for eqtl in eqtls:
            if eqtl[0]==eqtl_egene:
                for cell in cells:
                    cell = cell.strip()
                    to_file.append((eqtl_egene, eqtl[6], eqtl[3], row[5], cell))
                    print(eqtl_egene, eqtl[6], eqtl[3], row[5],cell)

    hic_out=open('hic_all_celltype.txt', 'w')
    writer = csv.writer(hic_out, delimiter='\t')
    writer.writerow(['Pair','Tissue','Pval','HiC'])
    writer.writerows(to_file)
    print(len(to_file))

if __name__=='__main__':
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-t', '--traits_dir',
                        default='../../results/trait_eqtls/',
                        help='Directory of trait directories containing'+\
                        'siginificant_eqtls.txt')
    parser.add_argument('-e', '--eqtlsDB',
                        default='../../results/spatial_eqtls.db',
                        help='Database of trait eQTLs and interactions.')
    parser.add_argument('-g', '--gwas_fp',
                        default='../../data/downloads/gwas_catalog_v1.0.1'+\
                        '-associations_e85_r2016-08-25.tsv',
                        help='Filepath to GWAS Catalog associations.')
    parser.add_argument('-o', '--output',
                        default='../../results/',
                        help='Filepath to write output.')
    args = parser.parse_args()
    #hic = 'interactions/interactions_summary.txt'
    gwas_genes = parse_gwas(args.gwas_fp)
    eqtls = parse_files(args.traits_dir, gwas_genes)
    #parse_hic(hic, eqtls)
