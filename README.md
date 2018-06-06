# Multimorbidity Atlas
To investigate trait or phenotype multimorbidity by integrating chromatin interaction and expression quantitative trait loci data.

### Data
- [GWAS Catalog](https://www.ebi.ac.uk/gwas) all associations v1.0.1 from [here](https://www.ebi.ac.uk/gwas/docs/file-downloads)
- [Rao et al (2014)](https://dx.doi/10.1016/j.cell.2014.11.021)'s Hi-C libraries of GM12878, IMR90, HMEC, NHEK, K562, HUVEC,  and KBM7 cell lines from the [GEO repository](https://ftp.ncbi.nih.gov/geo/samples/GSM1551nnn)
- [GTEx](https://www.gtexportal.org/) multi-tissue eQTLs data analysis v4
- [OMIM's](https://www.omim.org/) genemap2 data


## Getting Started 
- To begin, clone this repository with:
  ```
  git clone https://github.com/Genome3d/multimorbidity-atlas.git
  ```

- Install the [CoDeS3D][https://github.com/alcamerone/codes3d] pipeline:
  ```
  cd scripts/python/
  git clone https://github.com/alcamerone/codes3d.git
  ```


## Method overview

![Pipeline flow](docs/codeflow.png?raw=true "pipeline flow")


1. Obtaining SNPs from GWAS Catalog associations file from the command line:
   ```
   cd scripts/python/
   ./extract_gwas_snps.py
   ```
2. Identification of spatial eQTLs from GWAS SNPs:
   ```
   cd ../bash/
   ./batch_codes3d.sh
   ```
3. Saving eQTLs to database:
   ```
   cd ../python/
   ./init_db.py
   ```
4. Getting eQTLs for complex traits:
   ```
   ./get_trait_eqtls.py
   ```
5. Obtaining of tissue specific spatial eQTL-eGene interactions:
   ```
   cd ../bash/
   ./produce_summaries.sh
   ```
6. Construct matrices of shared eQTLs, eGenes ratios among complex traits. Plus control analysis:
   ```
   cd ../python/
   ./control_analysis.py
   ```
7. Cluster complex traits by shared eGenes and eQTLs
   ```
   cd ../R/
   ./convex_biclustering.R
   ```
8. Find traits that share eGene, using trait, SNPs, or genes as input:
   ```
   cd ../python/
   ./get_cluster.py
   ```