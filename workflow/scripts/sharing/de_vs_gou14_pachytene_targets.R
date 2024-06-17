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

# gou 2014
upInshMiwi<- readxl::read_xlsx("~/Downloads/41422_2014_BFcr201441_MOESM36_ESM.xlsx",sheet = 1) |>
  mutate(miwi.target=T)


dat |>
  dplyr::select(celltype,feature,gene_name,logFC) |>
  mutate(upInshMiwi = gene_name %in% upInshMiwi$GeneSymbol) |>
  ggplot(aes(upInshMiwi,logFC)) +
  geom_boxplot() +
  facet_wrap(~celltype,nrow=1) +
  ggpubr::stat_compare_means()



ranks <- split(dat,dat$celltype) |> map(~{dplyr::select(.x,gene_name,logFC)}) |> map(deframe)

gs <- list(gou14=upInshMiwi$GeneSymbol)
gsea_res <- ranks |>
  map(~fgsea(gs,stats=.x,nproc=1,eps=0))

gsea_res |>
  map_df(as_tibble,.id="label") |>
  mutate(padj = p.adjust(pval,method="BH"))  |>
  dplyr::select(label,pathway,pval,padj,NES) |>
  gt::gt()
