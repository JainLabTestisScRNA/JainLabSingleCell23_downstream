library(tidyverse)
library(scater)
library(scran)
library(scuttle)
library(TSCAN)

fl <- ifelse(exists("snakemake"),snakemake@input$sce,"results/germ_cells/adult.sce.integrated.clustered.celltypes.germ_cell.preprocessed.cell_cycle.rds")
sce <- read_rds(fl)

entropy <- perCellEntropy(sce,assay.type="logcounts")

sce$entropy <- entropy

# put x axis into the correct order
labs <- sce$label |> unique() 
level_order <- labs[as.integer(str_extract(labs,"\\d+")) |> order()]

sce$label <- sce$label |> fct_relevel(level_order)

g <- plotColData(sce, x = "genotype", y="entropy",other_fields = "label") +
  stat_summary() +
  facet_wrap(~label,) +
  theme(axis.text.x = element_text(angle=45, hjust=1))

g_entropy
