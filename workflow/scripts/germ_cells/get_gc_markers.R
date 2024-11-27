Sys.setenv(R_PROFILE=".Rprofile")
source(Sys.getenv("R_PROFILE"))

library(tidyverse)
library(scater)
library(scran)
library(scuttle)

sce <- read_rds(ifelse(exists("snakemake"),snakemake@input$sce,"results/germ_cells/adult.sce.germ_cell.wt.subclustered.reintegrated.rds"))

m.out <- scoreMarkers(sce, colLabels(sce), block=sce$batch)

m.out <- m.out |>
  map_df(as_tibble,rownames="feature",.id="label")

m.out <- rowData(sce)[,"gene_name",drop=F] |>
  as_tibble(rownames="feature") |>
  left_join(m.out,y=_,by="feature") |>
  dplyr::relocate(gene_name,.after="feature")

write_tsv(m.out,snakemake@output$tsv)
