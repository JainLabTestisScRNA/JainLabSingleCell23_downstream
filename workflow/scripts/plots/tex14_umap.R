library(tidyverse)
library(scater)
library(scran)
library(scuttle)

fl <- ifelse(exists("snakemake"),snakemake@input$sce,"results/integration/adult.sce.integrated.clustered.celltypes.rds")
sce <- read_rds(fl)

sce

sce <- sce[,order(logcounts(sce)["ENSMUSG00000010342",])]

g <- plotReducedDim(sce,"UMAP",ncomponents = 2,colour_by = "Tex14",swap_rownames = "gene_name",text_by = "celltype")



ggsave(snakemake@output$pdf,g)

g$data |> write_tsv(snakemake@output$tsv)
