library(tidyverse)
library(scater)
library(scran)
library(scuttle)
library(paletteer)

theme_set(theme_classic())

fl <- ifelse(exists("snakemake"),snakemake@input$sce,"results/germ_cells/adult.sce.integrated.clustered.celltypes.germ_cell.reprocessed.rds")
sce <- read_rds(fl)

x <- makePerCellDF(sce) |>
  as_tibble(rownames="cell") |>
  dplyr::select(cell,label,genotype,batch) |>
  mutate(label = fct_reorder(label,as.integer(str_extract(label,"\\d+"))))


g <- ggplot(x,aes(genotype,fill=label)) +
  geom_bar(position="fill") +
  scale_fill_paletteer_d("Polychrome::alphabet")
  
ggsave(snakemake@output$pdf,g)

g$data |> write_tsv(snakemake@output$tsv)