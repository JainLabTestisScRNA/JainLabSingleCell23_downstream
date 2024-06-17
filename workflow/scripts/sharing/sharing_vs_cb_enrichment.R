Sys.setenv(R_PROFILE=".Rprofile")
source(Sys.getenv("R_PROFILE"))

library(tidyverse)
library(scater)
library(scran)
library(scuttle)

dat <- read_rds(ifelse(exists("snakemake"),
                       snakemake@input$rds,
                       "results/sharing/sharing_score.rds"))

# meikar
cb_fold_enr <- readxl::read_xlsx("data/meikar2018.xlsx",skip = 2,sheet = 6) |>
  dplyr::select(gene_name ="name",cb.fold.enr="Fold enrichment")


# note that the supp data here is not a omplete list - only ~900 genes
dat |>
  pivot_wider(names_from = genotype,values_from = prop.zero) |>
  mutate(score = log2(MUT/WT)) |>
  left_join(cb_fold_enr) |>
  drop_na() |>
  ggplot(aes(score,cb.fold.enr)) +
  geom_point() +
  ggdensity::geom_hdr_points() +
  facet_wrap(~label,scales = "free")
