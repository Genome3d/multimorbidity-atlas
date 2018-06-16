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

- To quickly find multimorbid traits without running the entire methods,
download and unzip the [spatial eQTL data for traits](https://doi.org/10.17608/k6.auckland.6459728) into the results directory (1.59 GB). Then run step 9 in Method overview section below.

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
4. Get eQTLs associated with complex traits:
   ```
   ./get_trait_eqtls.py
   ```
5. Obtain significant tissue-specific spatial eQTLs:
   ```
   cd ../bash/
   ./produce_summaries.sh
   ```
6. Get eQTL-eGene interactions in tissues
   ```
   cd ../python/
   ./get_eqtl_interactions.py
   ```
7. Construct matrices of shared eQTLs, eGenes ratios among complex traits. Plus control analysis:
   ```
   ./control_analysis.py
   ```
8. Cluster complex traits by shared eGenes and eQTLs
   ```
   cd ../R/
   ./convex_biclustering.R
   ```
9. Find traits that share eGenes. Run command with -h for help:
   ```
   cd ../python/
   ./get_cluster.py
   ```