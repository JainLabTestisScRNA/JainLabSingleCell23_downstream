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


cb_enr <- readxl::read_xlsx("data/meikar2018.xlsx",skip = 2) |>
  pull(`Gene Symbol`) |>
  str_extract("\\w+")

dat |>
  dplyr::select(celltype,feature,gene_name,logFC) |>
  mutate(grr = gene_name %in% cb_enr) |>
  ggplot(aes(grr,logFC)) +
  geom_boxplot() +
  facet_wrap(~celltype,scales = "free") +
  ggpubr::stat_compare_means()
