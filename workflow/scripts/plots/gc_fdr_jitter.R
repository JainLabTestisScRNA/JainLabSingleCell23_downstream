Sys.setenv(R_PROFILE=".Rprofile")
source(Sys.getenv("R_PROFILE"))

library(tidyverse)

theme_set(theme_classic())

df <- read_tsv(ifelse(exists("snakemake"),snakemake@input$tsv,"results/differential_expression/adult.tbl.germ_cell.reprocessed.de.tsv.gz"))

g <- df |>
  mutate(celltype= fct_reorder(celltype,as.integer(str_extract(celltype,"\\d+")))) |>
  ggplot(aes(celltype,-log10(FDR))) +
  geom_point(position = position_jitter(seed=2),size=0.2) +
  geom_violin(scale = "width") +
  theme(axis.text.x=element_text(angle=45,hjust=1))

ggsave(snakemake@output$pdf, g,height=4)
  
write_tsv(g$data,snakemake@output$tsv)

  