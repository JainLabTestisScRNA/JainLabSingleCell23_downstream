library(tidyverse)
library(scater)
library(scran)
library(scuttle)
library(rtracklayer)
library(DropletUtils)


# ------------------------------------------------------------------------------
# organize input files

# read the main dataset we'll use
# this is raw matrices in h5 format
raw_fl <- Sys.glob("data/cellranger/custom_adult/adult1/raw_feature_bc_matrix.h5")
raw_fl <-  snakemake@input$raw

filt_fl <- Sys.glob("data/cellranger/custom_adult/adult1/filtered_feature_bc_matrix.h5")
filt_fl <-  snakemake@input$filt

genes_fl <- "data/genes.gtf.gz"
genes_fl <- snakemake@input$genes
# ------------------------------------------------------------------------------
# import
names(raw_fl) <- str_extract(raw_fl, "adult\\d+")
names(filt_fl) <- str_extract(filt_fl, "adult\\d+")

sce <- read10xCounts(filt_fl, col.names = T)
raw <- read10xCounts(raw_fl, col.names = T)

# ------------------------------------------------------------------------------
# add feature data
gene.data <- rtracklayer::import(genes_fl)

# this gets host genes and our custom repeats
gene.data <- gene.data[gene.data$type =="gene" | gene.data$gene_biotype %in% "repeat_element"]

names(gene.data) <- gene.data$gene_id

stopifnot(anyDuplicated(names(gene.data)) == 0)

# this commented section was used to identify the cmo features as the offenders
#sce <- sce[intersect(rownames(sce),names(gene.data)),]
#gene.data <- gene.data[intersect(rownames(sce),names(gene.data)),]

stopifnot(length(gene.data) == nrow(sce))

stopifnot(length(gene.data) == nrow(raw))

rowRanges(sce) <- gene.data[rownames(sce)]
rowRanges(raw) <- gene.data[rownames(raw)]

# ------------------------------------------------------------------------------
# export

saveRDS(sce, snakemake@output$rds)
saveRDS(raw, snakemake@output$raw)
