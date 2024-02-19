library(tidyverse)
library(scater)
library(scran)
library(scuttle)
library(edgeR)


sce_fl <- "results/single-cell-preproc/integrated/adult.cellranger.sce.integrated.rds"
sce_fl <- snakemake@input$sce
sce <- read_rds(sce_fl)

summed <- aggregateAcrossCells(sce, 
                               id=colData(sce)[,c("batch","genotype","celltype")])

summed$genotype <- summed$genotype |> factor() |> relevel(ref = "WT")

summed.filt <- summed[,summed$ncells >= 5]

de.results <- pseudoBulkDGE(summed.filt, 
                            label=summed.filt$celltype,
                            design=~batch + genotype,
                            coef="genotypeMUT",
                            condition=summed.filt$genotype ,
                            include.intermediates=T
)


df <- de.results |>
  map_df(as_tibble,rownames="feature",.id="celltype") |>
  filter(!is.na(logFC))

write_tsv(df,snakemake@output$tsv)
write_rds(summed,snakemake@output$pseudobulk)
write_rds(de.results,snakemake@output$de_results)
