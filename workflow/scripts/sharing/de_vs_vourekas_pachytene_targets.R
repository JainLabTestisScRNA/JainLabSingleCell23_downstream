Sys.setenv(R_PROFILE=".Rprofile")
source(Sys.getenv("R_PROFILE"))

library(tidyverse)
library(scater)
library(scran)
library(scuttle)
library(fgsea)

sce <- read_rds(ifelse(exists("snakemake"),
                       snakemake@input$sce,
                       "results/germ_cells/adult.sce.germ_cell.both_genotypes.subclustered.reintegrated.rds"))

lkup <- as_tibble(rowData(sce)[c("gene_name","gene_id")]) |>
  dplyr::rename(feature="gene_id")


dat <- read_tsv(ifelse(exists("snakemake"),
                       snakemake@input$tsv,
                       "results/differential_expression/adult.tbl.germ_cell.reprocessed.de.tsv.gz")) |>
  left_join(lkup)

# vourekas et al
miwi_targets <- readxl::read_xlsx("~/Downloads/41594_2012_BFnsmb2347_MOESM29_ESM.xlsx")

dat |>
  dplyr::select(celltype,feature,gene_name,logFC) |>
  mutate(miwi.target = gene_name %in% miwi_targets$symbol) |>
  ggplot(aes(miwi.target,logFC)) +
  geom_boxplot() +
  facet_wrap(~celltype, nrow=1) +
  ggpubr::stat_compare_means()


# https://bioconductor.org/packages/release/bioc/vignettes/fgsea/inst/doc/geseca-tutorial.html

ranks <- split(dat,dat$celltype) |> map(~{dplyr::select(.x,gene_name,logFC)}) |> map(deframe)

gs <- list(vourekas=miwi_targets$symbol)
gsea_res <- ranks |>
  map(~fgsea(gs,stats=.x,nproc=1,eps=0))

gsea_res |>
  map_df(as_tibble,.id="label") |>
  mutate(padj = p.adjust(pval,method="BH"))  |>
  dplyr::select(label,pathway,pval,padj,NES) |>
  gt::gt()

