Sys.setenv(R_PROFILE=".Rprofile")
source(Sys.getenv("R_PROFILE"))

library(tidyverse)
library(scater)
library(scran)
library(scuttle)
library(clusterProfiler)
library(DOSE)
library(enrichR)
library(enrichplot)
library(org.Mm.eg.db)

sce <- read_rds(ifelse(exists("snakemake"),snakemake@input$sce,"results/germ_cells/adult.sce.germ_cell.both_genotypes.subclustered.reintegrated.rds"))

lkup <- tibble(feature=rownames(sce),gene_name = rowData(sce)$gene_name)

de <- read_tsv("results/differential_expression/adult.tbl.germ_cell.reprocessed.de.tsv.gz")

de <- de |> left_join(lkup,by="feature")

universes <- de |>
  filter(str_detect(feature,"ENSMUSG")) |>
  group_by(celltype) |>
  summarise(genes=list(feature)) |>
  deframe()

sig.up <- de |>
  filter(FDR<0.05 & logFC >0) |>
  filter(str_detect(feature,"ENSMUSG")) |>
  group_by(celltype) |>
  summarise(genes=list(feature)) |>
  deframe()

sig.dn <- de |>
  filter(FDR<0.05 & logFC < 0) |>
  filter(str_detect(feature,"ENSMUSG")) |>
  group_by(celltype) |>
  summarise(genes=list(feature)) |>
  deframe()

random <- de |>
  filter(str_detect(feature,"ENSMUSG")) |>
  group_by(celltype) |>
  sample_n(500) |>
  summarise(genes=list(feature)) |>
  deframe()

write_tsv(tibble(sig.up$`4/Spermatocyte`),"~/Downloads/4sperm.up.txt",col_names = F)
write_tsv(tibble(sig.dn$`4/Spermatocyte`),"~/Downloads/4sperm.dn.txt",col_names = F)
write_tsv(tibble(universes$`4/Spermatocyte`),"~/Downloads/4sperm.univ.txt",col_names = F)

write_tsv(tibble(sig.up$`5/RoundSpermatid`),"~/Downloads/5round.up.txt",col_names = F)
write_tsv(tibble(sig.dn$`5/RoundSpermatid`),"~/Downloads/5round.dn.txt",col_names = F)
write_tsv(tibble(universes$`5/RoundSpermatid`),"~/Downloads/5round.univ.txt",col_names = F)

# manually done with DAVID:
# get "enrichments", but without background, 
# so this is only a way of annotating some interesting genes by their go
# term annotations

lkup_list <- lkup |> deframe()

read_tsv("~/Downloads/4sperm_dn.david.no_background.gobp_direct.txt") |>
  mutate(Genes = map(Genes,~unlist(str_split(.x,", ")))) |>
  mutate(Genes = map(Genes, ~{lkup_list[.x]})) |>
  mutate(Genes = map_chr(Genes,~paste(.x,collapse=", "))) |>
  write_tsv("~/Downloads/4sperm_dn.david.no_background.gobp_direct.gene_symbols.txt")

read_tsv("~/Downloads/5round_dn.david.no_background.gobp_direct.txt") |>
  mutate(Genes = map(Genes,~unlist(str_split(.x,", ")))) |>
  mutate(Genes = map(Genes, ~{lkup_list[.x]})) |>
  mutate(Genes = map_chr(Genes,~paste(.x,collapse=", "))) |>
  write_tsv("~/Downloads/5round_dn.david.no_background.gobp_direct.gene_symbols.txt")
