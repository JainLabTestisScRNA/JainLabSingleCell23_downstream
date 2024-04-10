library(tidyverse)
library(scater)
library(scran)
library(scuttle)
library(edgeR)


sce_fl <- "results/integration/adult.sce.integrated.clustered.celltypes.rds"
sce_fl <- snakemake@input$sce
sce <- read_rds(sce_fl)

sce$seqrun <- sce$batch |> str_remove("_\\d$")

summed <- aggregateAcrossCells(sce, 
                               id=colData(sce)[,c("seqrun","batch","genotype","celltype")])

summed$genotype <- summed$genotype |> factor() |> relevel(ref = "WT")

summed.filt <- summed[,summed$ncells >= 5]

summed.filt <- summed.filt[rowSums(counts(summed.filt)) > 5*ncol(summed.filt),]

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
