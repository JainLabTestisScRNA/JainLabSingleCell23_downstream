sSys.setenv(R_PROFILE=".Rprofile")
source(Sys.getenv("R_PROFILE"))

library(tidyverse)
library(scater)
library(scran)
library(scuttle)

sce <- read_rds(ifelse(exists("snakemake"),snakemake@input$sce,
                       "results/germ_cells/adult.sce.germ_cell.both_genotypes.subclustered.reintegrated.rds"))


colData(sce)


plotColData(sce,x = "subsets_chr9_percent",y = "subsets_chrX_percent",other_fields = c("celltype","genotype","label")) +
  facet_wrap(genotype~label,nrow=2,scales="free")

colData(sce) |>
    as_tibble(rownames="cell") |>
    dplyr::select(cell,genotype,label,total,nexprs,matches("chr") & contains("percent")) |>
    pivot_longer(-c(cell,genotype,label,total,nexprs)) |>
    ggplot(aes(label,value,fill=genotype)) +
    geom_boxplot() +
    facet_wrap(~name,scales = "free_y")
    

colData(sce) |>
  as_tibble(rownames="cell") |>
  dplyr::select(cell,genotype,label,total,nexprs,matches("chr")) |>
  group_by(genotype,label) |>
  pivot_longer(-c(cell,genotype,label,total,nexprs)) |>
  group_by(genotype,label,name) |>
  arrange(genotype,label,name,-value) |>
  mutate(value = row_number()) |>
  ungroup() |>
  pivot_wider(names_from = name,values_from = value) |>
  ggplot(aes(subsets_chrY_sum,subsets_chrX_sum)) +
  ggdensity::geom_hdr_points(size=0.1) +
  facet_wrap(genotype~label,nrow=2,scales="free")