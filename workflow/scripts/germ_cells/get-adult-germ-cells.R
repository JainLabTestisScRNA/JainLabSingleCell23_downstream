library(tidyverse)
library(scater)
library(scran)
library(scuttle)

fl <- "results/integration/adult.sce.integrated.clustered.celltypes.rds"
fl <- snakemake@input$sce

germ_types <- c("Spermatogonia","Spermatocyte","RoundSpermatid","Elongating")

sce <- read_rds(fl)

sce2 <- sce[,sce$celltype %in% germ_types]

sce2.mut <- sce2[,sce2$genotype == "MUT"]
sce2.wt <- sce2[,sce2$genotype == "WT"]

write_rds(sce2,snakemake@output$rds)
write_rds(sce2.mut,snakemake@output$mut_rds)
write_rds(sce2.wt,snakemake@output$wt_rds)