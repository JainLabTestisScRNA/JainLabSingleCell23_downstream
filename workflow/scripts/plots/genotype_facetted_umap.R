Sys.setenv(R_PROFILE=".Rprofile")
source(Sys.getenv("R_PROFILE"))

library(tidyverse)
library(scater)
library(scran)
library(scuttle)

theme_set(theme_classic())

fl <- ifelse(exists("snakemake"),snakemake@input$sce,"results/integration/adult.sce.integrated.clustered.celltypes.rds")
sce <- read_rds(fl)

g <- plotReducedDim(sce,"UMAP",ncomponents = 2,colour_by = "celltype",swap_rownames = "gene_name",text_by = "celltype",other_fields = "genotype")

g <- g + facet_wrap(~genotype)

ggsave(snakemake@output$pdf,g)

g$data |> write_tsv(snakemake@output$tsv)