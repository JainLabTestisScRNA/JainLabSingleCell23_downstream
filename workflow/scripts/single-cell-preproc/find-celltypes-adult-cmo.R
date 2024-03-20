library(batchelor)
library(scran)
library(scuttle)
library(tidyverse)
library(DropletUtils)
library(scater)
library(ggdensity)
library(scDblFinder)
library(flexmix)
library(celda)
library(bluster)
library(monocle3)
library(garnett)
library(org.Mm.eg.db)

fl <- snakemake@input$sce

corrected <- read_rds(fl)

garnett_fl <- "data/green2018_garnett_adult_testis_classifier.txt"
garnett_fl <- snakemake@input$garnett_fl

# ------------------------------------------------------------------------------
# cluster naming 
# https://cole-trapnell-lab.github.io/garnett/docs_m3/#loading-your-data
# https://cole-trapnell-lab.github.io/garnett/docs_m3/#installing-garnett

mono_cds <- new_cell_data_set(logcounts(corrected),cell_metadata = colData(corrected),gene_metadata = rowData(corrected))



marker_check <- check_markers(mono_cds, garnett_fl,
                              db=org.Mm.eg.db,
                              cds_gene_id_type = "ENSEMBL",
                              marker_file_gene_id_type = "SYMBOL")

g_plot_markers <- plot_markers(marker_check)

set.seed(20)
threads <- 8
threads <- snakemake@threads
garnett_classifier <- train_cell_classifier(cds = mono_cds,cores = threads,
                                            marker_file = garnett_fl,
                                            db=org.Mm.eg.db,
                                            cds_gene_id_type = "ENSEMBL",
                                            marker_file_gene_id_type = "SYMBOL")

mono_cds <- classify_cells(mono_cds, garnett_classifier,
                           db = org.Mm.eg.db,
                           cluster_extend = TRUE,
                           cds_gene_id_type = "ENSEMBL")


stopifnot(all(rownames(colData(mono_cds)) == rownames(colData(corrected))))

corrected$celltype <- colData(mono_cds)$cluster_ext_type

# ------------------------------------------------------------------------------
# cluster naming by vote - where some cells in a cluster have garnett labels that
# disagree with the predominant cluster type
colData(corrected) <- colData(corrected) |>
  as_tibble(rownames="cell") |>
  filter(celltype!="Unknown") |>
  dplyr::count(label, celltype) |>
  filter(!is.na(celltype)) |>
  group_by(label) |>
  slice_max(n,n=1,with_ties = F) |>
  ungroup() |> 
  dplyr::select(label, celltype) |>
  left_join(x=as_tibble(colData(corrected),rownames="cell"),y=_, by="label", suffix=c('.garnett',"")) |>
  column_to_rownames("cell") |>
  DataFrame()


corrected$label2 <- paste(corrected$label,corrected$celltype,sep="/")


write_rds(corrected, snakemake@output$rds)
write_rds(garnett_classifier, snakemake@output$garnett_classifier)
write_rds(marker_check, snakemake@output$marker_check)