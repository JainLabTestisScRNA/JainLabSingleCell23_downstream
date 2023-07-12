library(tidyverse)
library(rtracklayer)

fa_fl <- "data/dfam_sequences.fasta"
fa_fl <- snakemake@input$fasta

fa <- import(fa_fl)

res_fl <- "results/tradeseq/juvenile_13d_wt_null.cellranger.germ_cell.deg.rds"

res_fl <- snakemake@input$deg

conditionRes <- read_rds(res_fl)

dtes <- conditionRes$feature[conditionRes$padj <= 0.1 & !str_detect(conditionRes$feature,"^ENSMUSG")]

dtes <- fa[dtes]

export(dtes, snakemake@output$fa, format="fasta")

