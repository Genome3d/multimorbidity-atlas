# Multimorbidity Atlas
To investigate trait or phenotype multimorbidity by integrating chromatin interaction and expression quantitative trait loci data.

### Data
- [GWAS Catalog](https://www.ebi.ac.uk/gwas) all associations v1.0.1 from [here](https://www.ebi.ac.uk/gwas/docs/file-downloads)
- [Rao et al (2014)'s](https://dx.doi/10.1016/j.cell.2014.11.021) Hi-C libraries of GM12878, IMR90, HMEC, NHEK, K562, HUVEC,  and KBM7 cell lines from the [GEO repository](https://ftp.ncbi.nih.gov/geo/samples/GSM1551nnn)
- [GTEx](https://www.gtexportal.org/) multi-tissue eQTLs data analysis v4
- [OMIM's](https://www.omim.org/) genemap2 data


##Getting Started 
- To begin, clone this repository with:

...```
...git clone https://github.com/Genome3d/multimorbidity-atlas.git
...```

- Install the [CoDeS3D][https://github.com/alcamerone/codes3d] pipeline:
...```
...cd scripts/python/
...git clone https://github.com/alcamerone/codes3d.git
...```

## Method overview

![Pipeline flow](docs/codeflow.png?raw=true "pipeline flow")


1. Extract SNPs from GWAS Catalog associations file from the command line:
...```
...cd scripts/python/
...python extract_gwas_snps.py
...```

