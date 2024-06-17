Sys.setenv(R_PROFILE=".Rprofile")
source(Sys.getenv("R_PROFILE"))

library(tidyverse)
library(scater)
library(scran)
library(scuttle)
library(fgsea)

dat <- read_rds(ifelse(exists("snakemake"),
                       snakemake@input$rds,
                       "results/sharing/sharing_score.rds"))


upInshMiwi<- readxl::read_xlsx("~/Downloads/41422_2014_BFcr201441_MOESM36_ESM.xlsx",sheet = 1)


dat |>
  pivot_wider(names_from = genotype,values_from = prop.zero) |>
  mutate(score = log2(MUT/WT)) |>
  mutate(miwi.target = gene_name %in% upInshMiwi$GeneSymbol) |>
  ggplot(aes(miwi.target,score)) +
  geom_boxplot() +
  facet_wrap(~label,scales = "free")



ranks <- split(dat,dat$celltype) |> map(~{dplyr::select(.x,gene_name,logFC)}) |> map(deframe)

gs <- list(goh15=pachy_pirna_targets$Gene)
gsea_res <- ranks |>
  map(~fgsea(gs,stats=.x,nproc=1,eps=0))

gsea_res |>
  map_df(as_tibble,.id="label") |>
  mutate(padj = p.adjust(pval,method="BH"))  |>
  dplyr::select(label,pathway,pval,padj,NES) |>
  gt::gt()