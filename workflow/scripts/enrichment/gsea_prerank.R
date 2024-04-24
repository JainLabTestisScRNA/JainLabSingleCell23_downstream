Sys.setenv(R_PROFILE=".Rprofile")
source(Sys.getenv("R_PROFILE"))

library(tidyverse)
library(scater)
library(scran)
library(scuttle)
library(fgsea)
library(GSEABase)

sce <- read_rds(ifelse(exists("snakemake"),snakemake@input$sce,"results/germ_cells/adult.sce.germ_cell.both_genotypes.subclustered.reintegrated.rds"))

lkup <- tibble(feature=rownames(sce),gene_name = rowData(sce)$gene_name)

de <- read_tsv("results/differential_expression/adult.tbl.germ_cell.reprocessed.de.tsv.gz")

de <- de |> left_join(lkup,by="feature")

ranks <- split(de,de$celltype) |> map(~{dplyr::select(.x,gene_name,logFC)}) |> map(deframe)

# https://bioconductor.org/packages/release/bioc/vignettes/fgsea/inst/doc/geseca-tutorial.html

gmt <- gmtPathways(ifelse(exists("snakemake"),snakemake@input$msigdb,"data/m5.all.v2023.2.Mm.symbols.gmt"))

gseaRes <- map(ranks, ~fgsea(gmt,stats = .x, minSize = 10, maxSize = 300,nproc=1))


write_rds(gseaRes, snakemake@output$rds)

