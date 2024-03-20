library(tidyverse)
library(scater)
library(scran)
library(scuttle)

fl <- "results/integration/adult.sce.integrated.clustered.celltypes.rds"
fl <- snakemake@input$sce

germ_types <- c("Spermatogonia","Spermatocyte","RoundSpermatid","Elongating")

sce <- read_rds(fl)

sce2 <- sce[,sce$celltype %in% germ_types]

write_rds(sce2,snakemake@output$rds)