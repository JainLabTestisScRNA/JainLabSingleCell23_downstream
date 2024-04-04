library(tidyverse)
library(scater)
library(scran)
library(scuttle)

theme_set(theme_classic())

fl <- ifelse(exists("snakemake"),snakemake@input$sce,"results/germ_cells/adult.sce.integrated.clustered.celltypes.germ_cell.preprocessed.cell_cycle.rds")
sce <- read_rds(fl)

x <- makePerCellDF(sce) |>
  as_tibble(rownames="cell") |>
  dplyr::select(cell,label,genotype,batch,cell_cycle)


g <- x |>
  ggplot(aes(label,fill=cell_cycle)) +
  geom_bar(position="fill") +
  ylab("proportion")

# put x axis into the correct order
labs <- sce$label |> unique() 
level_order <- labs[as.integer(str_extract(labs,"\\d+")) |> order()]
g <- g + scale_x_discrete(limits = level_order)

ggsave(snakemake@output$pdf,g)

g$data |> write_tsv(snakemake@output$tsv)