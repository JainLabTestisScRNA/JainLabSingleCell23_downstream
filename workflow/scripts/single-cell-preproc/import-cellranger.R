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
# naming and merging duplicate rows

# # rename with gene symbol so it's easier to subset and visualize for humans.
# names(sce) <- sce@rowRanges$gene_name
# 
# # this produces features with duplicated names (ie 1 named gene with multiple features) n=~38 total gene names
# duplicated_names <- rownames(sce)[duplicated(rownames(sce))] |> unique() |> set_names()
# 
# # the integer ixs of the duplicates
# duplicated_ix <- duplicated_names |>
#   map(~{which(rownames(sce) == .x)})
# 
# # sce with non duplicated gene names
# non_duplicated_sce <- sce[!rownames(sce) %in% duplicated_names,]
# 
# # generate the sce subset for duplicated gene names. we arbitrarily take the
# # first range for each duplicated gene.s
# merged_matrix <- sumCountsAcrossFeatures(sce, ids=duplicated_ix)
# merged_rngs <-  map(duplicated_ix, ~{rowRanges(sce)[.x[1]]}) |> GRangesList() |> unlist()
# names(merged_rngs) <- merged_rngs$gene_name
# 
# # expectation checks
# stopifnot(all(rownames(merged_matrix) == names(merged_rngs)))
# stopifnot(all(colnames(merged_matrix) == colnames(sce)))
# 
# # make new sce with the merged features - 1 row/rng for each duplicated gene name
# merged_sce <- SingleCellExperiment(assays=list(counts=merged_matrix),
#                                    rowRanges = merged_rngs, 
#                                    colData=colData(sce))
# 
# 
# sce <- rbind(non_duplicated_sce, merged_sce) |> sort()
# 
# 
# 
# stopifnot(anyDuplicated(rownames(sce2))==0)

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
