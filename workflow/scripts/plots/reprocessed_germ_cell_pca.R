Sys.setenv(R_PROFILE=".Rprofile")
source(Sys.getenv("R_PROFILE"))

library(tidyverse)
library(scater)
library(scran)
library(scuttle)

theme_set(theme_classic())

fl <- ifelse(exists("snakemake"),snakemake@input$sce,"results/germ_cells/adult.sce.germ_cell.both_genotypes.subclustered.reintegrated.rds")
sce <- read_rds(fl)


g12 <- plotReducedDim(sce,"corrected",
                     ncomponents = c(1,2),
                     colour_by = "label",swap_rownames = "gene_name",text_by = "label") +
  coord_fixed()

g13 <- plotReducedDim(sce,"corrected",
                      ncomponents = c(1,3),
                      colour_by = "label",swap_rownames = "gene_name",text_by = "label") +
  coord_fixed()


g12_gt <- plotReducedDim(sce,"corrected",other_fields = "genotype",
                      ncomponents = c(1,2),
                      colour_by = "label",swap_rownames = "gene_name",text_by = "label") +
  facet_wrap(~genotype) +
  coord_fixed() 

g13_gt <- plotReducedDim(sce,"corrected", other_fields = "genotype",
                      ncomponents = c(1,3),
                      colour_by = "label",swap_rownames = "gene_name",text_by = "label") +
  facet_wrap(~genotype) +
  coord_fixed()

pdf(snakemake@output$pdf)
g12

g13

g12_gt

g13_gt

dev.off()


sce |> makePerCellDF() |>
  as_tibble(rownames="cell") |>
  dplyr::select(cell,genotype,batch,label,corrected.1,corrected.2,corrected.3) |>
  write_tsv(snakemake@output$tsv)
