Sys.setenv(R_PROFILE=".Rprofile")
source(Sys.getenv("R_PROFILE"))

library(tidyverse)
library(scater)
library(scran)
library(scuttle)

sce <- read_rds(ifelse(exists("snakemake"),
                       snakemake@input$sce,
                       "results/germ_cells/adult.sce.germ_cell.both_genotypes.subclustered.reintegrated.rds"))

lkup <- as_tibble(rowData(sce)[c("gene_name","gene_id")]) |>
  dplyr::rename(feature="gene_id")


dat <- read_tsv(ifelse(exists("snakemake"),
                       snakemake@input$tsv,
                       "results/differential_expression/adult.tbl.germ_cell.reprocessed.de.tsv.gz")) |>
  left_join(lkup)


grr <- readxl::read_xlsx("~/Downloads/41586_2018_BFnature25964_MOESM3_ESM.xlsx",col_names = "grr")$grr

dat |>
  dplyr::select(celltype,feature,gene_name,logFC) |>
  mutate(grr = gene_name %in% grr) |>
  ggplot(aes(grr,logFC)) +
  geom_boxplot() +
  facet_wrap(~celltype,scales = "free") +
  ggpubr::stat_compare_means()
