# Cancer-CRISPR
This repository contains the Nextflow pipeline for Cancer-CRISPR tool. This pipeline maps the sgRNA sequences provided in the library.fa file to hg38, retrieves mapped chromosome, strand as well as start and end positions. Further it checks whether the gene a particular sgRNA mapped to matches the provided gene name for the sgRNA. Finally, it downloads gene expression matrix from 2 TCGA-BRCA datasets and provides the TPM expression values obtained from STAR - Counts workflow. The mapping table (reads_mapped_table.tsv) and the TCGA expression data table (tcga_data_table.tsv) are created inside the results directory when the pipeline is run.

This pipeline uses a prebult index of hg38 which should be provided by the user inside data/index folder.

**Workflow:**
<p align="left">
  <img src="https://github.com/erkutilaslan/cancer_crispr/blob/main/dag.png"></div>
</p>


# Dependencies
This workflow is written in DSL1 so it is compatible with Nextflow version 22.10.x or earlier.

Required software to run the pipeline is provided in env_explicit.yml file. 

Conda environment using this file can be built with the command:
```bash
conda env create -f env_explicit.yml
```

**R packages:**

This workflow uses TCGAbiolinks package to retrieve gene expression data of TCGA and biomaRt to map ENSEMBL gene IDs to HGNC.
To install these packages run the following code in R environment:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("dplyr","TCGAbiolinks","biomaRt"))
```

# Usage
- Put pre-built index files inside the data/index/ directory.

- Run the pipeline:
```bash
nextflow run cancer_crispr.nf
```
  

