#!/usr/bin/env Rscript

# load package
library(TCGAbiolinks)
library(SummarizedExperiment)
library(biomaRt)
library(stringr)
library(dplyr)

# setup variable to control command line args
args <- commandArgs(TRUE)

# desired samples to be downloaded
listSamples <- c(
    "TCGA-A7-A13D-01A-13R-A12P-07", "TCGA-E9-A1RH-11A-34R-A169-07"
)


# query parameters
query <- GDCquery(
    project = "TCGA-BRCA",
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    experimental.strategy = "RNA-Seq",
    barcode = listSamples,
    workflow.type = "STAR - Counts",
    access = "open"
)

# download query
GDCdownload(query)

# prepare expression matrix with geneID in the rows and samples (barcode) in the columns
BRCA.Rnaseq.SE <- GDCprepare(query)
BRCAMatrix <- assay(BRCA.Rnaseq.SE, "tpm_unstrand")

# convert matrix into dataframe
brca_df <- data.frame(ensembl_id = row.names(BRCAMatrix), BRCAMatrix)

# map ensembl gene ids to hgnc
ensembl <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

# remove version numbers from ensembl_ids for better mapping
brca_df$ensembl_id <- str_replace(brca_df$ensembl_id,
                        pattern = ".[0-9]+$",
                        replacement = "")

#prepare ensembl id list for mapping
ensembl_ids <- brca_df$ensembl_id

ensembl_to_hgnc <- getBM(attributes = c('ensembl_gene_id',
                                        'hgnc_symbol'),
                         filters = 'ensembl_gene_id',
                         values = ensembl_ids,
                         mart = ensembl)

# add HGNC by join
brca_df <- left_join(brca_df, ensembl_to_hgnc,
                     by = c("ensembl_id" = "ensembl_gene_id"),
                     copy = FALSE)


# import the data from previous nextflow process
# i need to change this step when I convert this to Rscript
sgrna_data <- read.table(args[1], header = TRUE)
#sgrna_data <- read.table("../results/reads_mapped_table.tsv", header = TRUE)


# filter genes only present in sgrna_map_data
brca_df <- filter(brca_df, hgnc_symbol %in% sgrna_data$Target_gene)

# reorder brca_df
brca_df <- brca_df[, c(4,2,3)]

# remove NAs
brca_df <- na.omit(brca_df)

#output brca_df
write.table(
            brca_df,
            "tcga_data_table.tsv",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE
)


