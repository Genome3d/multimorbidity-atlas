#!/usr/bin/env python
import os 
import sqlite3
import csv

def createDB(db_fp):
    """Creates database to store all eQTLs in eqtls.txt files from 
    ../bash/batch_codes3d.sh
    Note that the initial run in GTEx V4 does not include eQTL effect sizes.

    Args:
        db_fp: Filepath of the SQLite database to be created.

    Output:
       spatial_eqtls.db: A db with an eqtls table.
    """
    conn = sqlite3.connect(os.path.join(db_fp,'spatial_eqtls.db'))
    conn.text_factory = str
    cur = conn.cursor()
    cur.execute("DROP TABLE IF EXISTS eqtls")
    cur.execute(""" 
        CREATE TABLE eqtls 
        (snp TEXT, snp_chr INTEGER, snp_locus INTEGER, 
            gene TEXT, gene_chr INTEGER, gene_start INTEGER, gene_end INTEGER,
            cell_line TEXT, interaction TEXT, 'snp-gene_distance' INTEGER, 
            tissue TEXT, pvalue REAL)
        """)
    conn.close()

def popDB(db_fp, results_dir_fp):
    """Creates database to store all eQTLs in eqtls.txt files from 
    batch_codes3d.sh run.
    Note that the initial run in GTEx V4 does not include eQTL effect sizes.

    Args:
        db_fp: Filepath of the SQLite database to be created.
        results_dir_fp: Filepath of ../bash/batch_codes3d.sh run

    Output:
       spatial_eqtls.db: with a populated eqtls table.
    """
    
    conn = sqlite3.connect(os.path.join(db_fp,'./spatial_eqtls.db'))
    conn.text_factory = str
    cur = conn.cursor()
    dirs = [os.path.join(results_dir_fp, d) for d in os.listdir(results_dir_fp)]
    for d in dirs:
        dfile = open(d + '/eqtls.txt', 'r')
        ddata = csv.reader(dfile, delimiter = '\t')
        next(ddata, None)
        for row in ddata:
            cur.executemany("""INSERT INTO eqtls  VALUES 
            (?,?,?,?,?,?,?,?,?,?,?,?);""", (row,))
    conn.commit()
    conn.close()



if __name__ == "__main__":
    results_dir_fp = '../../results/batch_results/'
    db_fp = '../../results/'
    createDB(db_fp)
    popDB(db_fp, results_dir_fp)

