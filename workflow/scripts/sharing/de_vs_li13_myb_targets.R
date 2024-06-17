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

# li et al 2013
# https://www.sciencedirect.com/science/article/pii/S109727651300172X?via%3Dihub#app3
myb_targets <- readxl::read_xlsx("data/li2013_table3_1-s2.0-S109727651300172X-mmc3.xlsx",skip = 1)

de_myb_targets <-  filter(myb_targets,`q-value` < 0.05) |> pull(Gene)
direct_myb_targets <- filter(myb_targets,`Distance to summit of nearest A-Myb ChIP-seq peak (bp)` < 10000) |> pull(Gene)


dat |>
  dplyr::select(celltype,feature,gene_name,logFC) |>
  mutate(myb.target = gene_name %in%direct_myb_targets) |>
  ggplot(aes(myb.target,logFC)) +
  geom_boxplot() +
  facet_wrap(~celltype) +
  ggpubr::stat_compare_means()

ranks <- split(dat,dat$celltype) |> map(~{dplyr::select(.x,gene_name,logFC)}) |> map(deframe)

gs <- list(li13=direct_myb_targets)
gsea_res <- ranks |>
  map(~fgsea(gs,stats=.x,nproc=1,eps=0))

gsea_res |>
  map_df(as_tibble,.id="label") |>
  mutate(padj = p.adjust(pval,method="BH"))  |>
  dplyr::select(label,pathway,pval,padj,NES) |>
  gt::gt()
