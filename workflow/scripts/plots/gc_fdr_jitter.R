library(tidyverse)

theme_set(theme_classic())

df <- read_tsv(ifelse(exists("snakemake"),snakemake@input$tsv,"results/differential_expression/adult.tbl.germ_cell.reprocessed.de.tsv.gz"))

g <- df |>
  mutate(celltype= fct_reorder(celltype,as.integer(str_extract(celltype,"\\d+")))) |>
  ggplot(aes(celltype,-log10(FDR))) +
  geom_point(position = position_jitter(seed=2)) +
  geom_violin(scale = "width")

ggsave(snakemake@output$pdf, g)
  
write_tsv(g$data,snakemake@output$tsv)

  