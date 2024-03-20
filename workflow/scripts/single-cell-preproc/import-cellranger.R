library(tidyverse)
library(scater)
library(scran)
library(scuttle)
library(rtracklayer)
library(DropletUtils)


# ------------------------------------------------------------------------------
# organize input files

# read the main dataset we'll use
# this is two filtered matrices in h5 format
# for consistency with the unfiltered matrix, I'll remove the 1_/2_ prefix
fl <- "data/cellranger/custom/raw_feature_bc_matrix/"
fl <-  snakemake@input$mat

genes_fl <- "data/genes.gtf.gz"
genes_fl <- snakemake@input$genes

# sample calls from cellranger multi pipeline
cmo_fl <- "data/cellranger/custom/assignment_confidence_table.csv"
cmo_fl <- snakemake@input$sample_assignments

# ------------------------------------------------------------------------------
# import
sce <- read10xCounts(fl, col.names = T)

# this can show up when importing multiple matrices at once
# keep for consistency, or if i repurpose this script
colnames(sce) <- str_remove(colnames(sce),"^\\d_")

# ------------------------------------------------------------------------------
# add feature data
gene.data <- rtracklayer::import(genes_fl)

# remove cmo counts, bc we rely on cellrangers sample
# assignments and don't want these for downstream stuff.
sce <- sce[!str_detect(rownames(sce),"CMO"),]

# this gets host genes and our custom repeats
gene.data <- gene.data[gene.data$type =="gene" | gene.data$gene_biotype %in% "repeat_element"]

names(gene.data) <- gene.data$gene_id

stopifnot(anyDuplicated(names(gene.data)) == 0)

# this commented section was used to identify the cmo features as the offenders
#sce <- sce[intersect(rownames(sce),names(gene.data)),]
#gene.data <- gene.data[intersect(rownames(sce),names(gene.data)),]

stopifnot(length(gene.data) == nrow(sce))

rowRanges(sce) <- gene.data[rownames(sce)]

# ------------------------------------------------------------------------------
# add sample data

cmo <- read_csv(cmo_fl)

newcoldata <- colData(sce) |>
  as_tibble(rownames="cell") |>
  left_join(cmo, by="Barcode", relationship="one-to-one") |>
  column_to_rownames("cell") |>
  DataFrame()

stopifnot(rownames(newcoldata) == rownames(colData(sce)))

colData(sce) <- newcoldata

sce$Sample <- sce$Assignment

# ------------------------------------------------------------------------------
# export

sce_raw <- sce

sce <- sce[,!is.na(sce$Sample) & !sce$Sample %in% c("Blank","Multiplet","Unassigned")]

saveRDS(sce, snakemake@output$rds)
saveRDS(sce_raw, snakemake@output$raw)
