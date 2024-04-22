library(tidyverse)
library(scater)
library(scran)
library(scuttle)
library(edgeR)


sce_fl <- "results/germ_cells/adult.sce.germ_cell.both_genotypes.subclustered.reintegrated.rds"
sce_fl <- snakemake@input$sce
sce <- read_rds(sce_fl)

sce$seqrun <- sce$batch |> str_remove("_\\d$")

summed <- aggregateAcrossCells(sce, 
                               id=colData(sce)[,c("seqrun","batch","genotype","label")])

summed$genotype <- summed$genotype |> factor() |> relevel(ref = "WT")

summed.filt <- summed[,summed$ncells >= 5]

summed.filt <- summed.filt[rowSums(counts(summed.filt)) > 5*ncol(summed.filt),]

de.results <- pseudoBulkDGE(summed.filt, 
                            label=summed.filt$label,
                            design=~batch + genotype,
                            coef="genotypeMUT",
                            method="edgeR",
                            condition=summed.filt$genotype ,
                            include.intermediates=T
)


df <- de.results |>
  map_df(as_tibble,rownames="feature",.id="celltype") |>
  #mutate(FDR = adj.P.Val) |>
  filter(!is.na(logFC))

write_tsv(df,snakemake@output$tsv)
write_rds(summed,snakemake@output$pseudobulk)
write_rds(de.results,snakemake@output$de_results)
