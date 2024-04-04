library(tidyverse)
library(scater)
library(scran)
library(scuttle)

theme_set(theme_classic())

fl <- ifelse(exists("snakemake"),snakemake@input$sce,"results/integration/adult.sce.integrated.clustered.celltypes.rds")

sce <- read_rds(fl)

g <- plotExpression(sce,x = "celltype",
                    features = c("Tex14"),
                    swap_rownames = "gene_name",scales = "free") +
  theme(axis.text.x=element_text(angle=45,hjust=1))


ggsave(snakemake@output$pdf,g)

g$data |> write_tsv(snakemake@output$tsv)
