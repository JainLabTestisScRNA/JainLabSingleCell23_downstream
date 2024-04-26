Sys.setenv(R_PROFILE=".Rprofile")
source(Sys.getenv("R_PROFILE"))

library(tidyverse)
library(scater)
library(scran)
library(scuttle)

theme_set(theme_classic())

fl <- ifelse(exists("snakemake"),snakemake@input$sce,"results/germ_cells/adult.sce.integrated.clustered.celltypes.germ_cell.reprocessed.rds")
sce <- read_rds(fl)

oi <- list(Sycp3="ENSMUSG00000020059",Stra8="ENSMUSG00000029848",Dmrt1="ENSMUSG00000024837")



g1 <- plotReducedDim(sce[,order(logcounts(sce)[oi$Sycp3,])],"corrected",
               ncomponents = c(1,3),
               colour_by = "Sycp3",swap_rownames = "gene_name",text_by = "label")

g2 <- plotReducedDim(sce[,order(logcounts(sce)[oi$Stra8,])],"corrected",
               ncomponents = c(1,3),
               colour_by = "Stra8",swap_rownames = "gene_name",text_by = "label")

g3 <- plotReducedDim(sce[,order(logcounts(sce)[oi$Dmrt1,])],"corrected",
               ncomponents = c(1,3),
               colour_by = "Dmrt1",swap_rownames = "gene_name",text_by = "label")

gl <- list(Sycp3=g1,Stra8=g2,Dmrt1=g3)

pdf(snakemake@output$pdf)
gl
dev.off()

gl |>
  map(`$`,"data") |>
  map_df(~as_tibble(.x,rownames="cell"),.id="gene") |>
  write_tsv(snakemake@output$tsv)
