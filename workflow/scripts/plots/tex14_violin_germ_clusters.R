library(tidyverse)
library(scater)
library(scran)
library(scuttle)

fl <- ifelse(exists("snakemake"),snakemake@input$sce,"results/germ_cells/adult.sce.integrated.clustered.celltypes.germ_cell.reprocessed.rds")

sce <- read_rds(fl)

g <- plotExpression(sce,x = "label",
                    features = c("Tex14"),
                    swap_rownames = "gene_name",scales = "free") +
  theme(axis.text.x=element_text(angle=45,hjust=1))

# put x axis into the correct order
labs <- sce$label |> unique() 
level_order <- labs[as.integer(str_extract(labs,"\\d+")) |> order()]
g <- g + scale_x_discrete(limits = level_order)

ggsave(snakemake@output$pdf,g)

g$data |> write_tsv(snakemake@output$tsv)
