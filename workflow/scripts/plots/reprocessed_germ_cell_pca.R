library(tidyverse)
library(scater)
library(scran)
library(scuttle)

theme_set(theme_classic())

fl <- ifelse(exists("snakemake"),snakemake@input$sce,"results/germ_cells/adult.sce.integrated.clustered.celltypes.germ_cell.reprocessed.rds")
sce <- read_rds(fl)


g12 <- plotReducedDim(sce,"corrected",
                     ncomponents = c(1,2),
                     colour_by = "label",swap_rownames = "gene_name",text_by = "label")

g13 <- plotReducedDim(sce,"corrected",
                      ncomponents = c(1,3),
                      colour_by = "label",swap_rownames = "gene_name",text_by = "label")

ggsave(snakemake@output$pc1_pc2,g12)
ggsave(snakemake@output$pc1_pc3,g13)


sce |> makePerCellDF() |>
  as_tibble(rownames="cell") |>
  dplyr::select(cell,genotype,batch,label,corrected.1,corrected.2,corrected.3) |>
  write_tsv(snakemake@output$tsv)
