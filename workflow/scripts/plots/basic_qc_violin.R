Sys.setenv(R_PROFILE=".Rprofile")
source(Sys.getenv("R_PROFILE"))

library(tidyverse)
library(scater)
library(scran)
library(scuttle)

x_by <- "label"
x_by <- snakemake@params$x

sce <- read_rds(ifelse(exists("snakemake"),snakemake@input$sce,
                       "results/germ_cells/adult.sce.germ_cell.both_genotypes.subclustered.reintegrated.rds"))

sce$celltype <- colData(sce)[,x_by]

g <- makePerCellDF(sce) |>
  as_tibble(rownames="cell") |>
  #mutate(across(contains("percent"),~{replace_na(na_if(.x,0),-Inf)})) |>
  dplyr::select(cell,celltype,genotype,sum,detected,subsets_MT_percent,subsets_TE_percent) |>
  pivot_longer(-c(cell,celltype,genotype)) |>
  mutate(celltype = fct_relevel(celltype,
                                c("Sertoli","Macrophage","Endothelial","InnateLymphoid","Spermatogonia","Spermatocyte","RoundSpermatid","Elongating"))) |>
  mutate(celltype=fct_rev(celltype)) |>
  mutate(name=fct_relevel(name,c("sum","detected"))) |>
  ggplot(aes(value,celltype)) +
  geom_violin(scale = "width") +
  stat_summary(geom = "point",color="red",fun = "mean") +
  facet_grid(genotype~name,scales="free",switch = "x") +
  scale_x_log10(oob=scales::oob_squish_infinite) +
  xlab("")

ggsave(snakemake@output$pdf, g)

write_tsv(g$data,snakemake@output$tsv)

