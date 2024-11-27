Sys.setenv(R_PROFILE=".Rprofile")
source(Sys.getenv("R_PROFILE"))

library(tidyverse)
library(scater)
library(scran)
library(scuttle)

dat <- read_rds(ifelse(exists("snakemake"),
                       snakemake@input$rds,
                       "results/sharing/sharing_score.rds"))


# goh et al
pachy_pirna_targets <- readxl::read_xlsx("~/Downloads/TableS7_628_low-confidence_mRNA_Targets.xlsx")

dat |>
  pivot_wider(names_from = genotype,values_from = prop.zero) |>
  mutate(score = log2(MUT/WT)) |>
  mutate(pachy = gene_name %in% pachy_pirna_targets$Gene) |>
  drop_na() |>
  ggplot(aes(pachy,score)) +
  geom_boxplot() +
  facet_wrap(~label,scales = "free")
