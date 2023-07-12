library(tidyverse)
library(scater)
library(scran)
library(scuttle)
library(ggpubr)
library(ggdensity)
library(patchwork)
library(paletteer)
library(dtwclust)
library(tradeSeq)
library(pheatmap)
library(TSCAN)
library(slingshot)

sce_fl <- "results/tradeseq/juvenile_13d_wt_null.cellranger.sce.germ_cell.tradeseq.rds"
sce_fl <- snakemake@input$sce

sce <- read_rds(sce_fl)

conditionRes <- conditionTest(sce)

id_lkup <- rowData(sce)[,c("gene_id", "gene_name")] |>
  as_tibble() |>
  dplyr::rename(feature="gene_id")

conditionRes <- conditionRes |> 
  as_tibble(rownames="feature") |> drop_na() |>
  mutate(padj = p.adjust(pvalue,method="BH")) |>
  arrange(-abs(waldStat)) |>
  arrange(pvalue) |>
  left_join(id_lkup, by="feature")

write_rds(conditionRes, file = snakemake@output$deg_rds)
