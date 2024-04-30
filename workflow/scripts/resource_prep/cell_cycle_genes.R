Sys.setenv(R_PROFILE=".Rprofile")
source(Sys.getenv("R_PROFILE"))

library(tidyverse)
library(jsonlite)
library(readxl)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(Orthology.eg.db)
library(msigdb)
library(ExperimentHub)
library(GSEABase)

# ------------------------------------------------------------------------------
# macosko et al cell cyce signature
# ------------------------------------------------------------------------------

genes <- read_xlsx(ifelse(exists("snakemake"),snakemake@input$xlsx,"data/macosko2015_supp2.xlsx"),sheet = 2) |>
  pivot_longer(everything()) |>
  arrange(name) |>
  drop_na()

# credit: https://support.bioconductor.org/p/9153600/
hg <- mapIds(org.Hs.eg.db, genes$value, "ENTREZID", "SYMBOL")

mapped <- select(Orthology.eg.db, hg, "Mus_musculus", "Homo_sapiens")

musymb <- select(org.Mm.eg.db, as.character(mapped[,2]), "SYMBOL","ENTREZID") |> as_tibble() |> distinct()
muid <- select(org.Mm.eg.db, as.character(mapped[,2]), "ENSEMBL","ENTREZID") |> as_tibble() |> distinct()

conv <- as_tibble(mapped) |>
  mutate(across(everything(), as.character)) |>
  distinct() |>
  left_join(enframe(hg,name = "SYMBOL",value = "ENTREZID"),by=c(Homo_sapiens="ENTREZID")) |>
  left_join(musymb,by=c(Mus_musculus="ENTREZID"),suffix=c(".hs",'.mm')) |>
  left_join(muid,by=c(Mus_musculus="ENTREZID")) |>
  distinct()
  
df <- left_join(genes,conv,by=c("value"="SYMBOL.hs"))

write_tsv(df,snakemake@output$tsv)