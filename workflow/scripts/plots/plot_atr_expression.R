Sys.setenv(R_PROFILE=".Rprofile")
source(Sys.getenv("R_PROFILE"))

library(tidyverse)
library(scater)
library(scran)
library(scuttle)

theme_set(theme_classic())

read_tsv("results/plots/mut_vs_wt_volcanos_and_ma/germ_cell_mut_vs_wt_volcano.tsv.gz") -> de
de |> filter(gene_name == "Atr")

fl <- ifelse(exists("snakemake"),snakemake@input$sce,"results/germ_cells/adult.sce.germ_cell.both_genotypes.subclustered.reintegrated.rds")
sce <- read_rds(fl)

g <- plotExpression(sce,x = "label",
               features = c("Atr"),
               swap_rownames = "gene_name",scales = "free",other_fields = "genotype") +
  theme(axis.text.x=element_text(angle=45,hjust=1))

g+ facet_wrap(~genotype)


ggsave(snakemake@output$pdf,g)

g$data |> write_tsv(snakemake@output$tsv)
