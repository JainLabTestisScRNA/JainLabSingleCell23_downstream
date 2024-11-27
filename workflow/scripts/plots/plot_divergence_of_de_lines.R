library(tidyverse)

de <- read_tsv("results/differential_expression/adult.tbl.germ_cell.reprocessed.de.tsv.gz")
div <- read_tsv("~/work/mouse_te_ages/mm10.div.clean.tsv")
de |>
  filter(celltype == "4/Spermatocyte") |>
  filter(str_detect(feature,"L1|L2|Lx")) |>
  left_join(div,by=c(feature="Repeat")) |>
  mutate(is.de = FDR<0.05) |>
  filter(Kimura < 0)
  ggplot(aes(is.de,Kimura)) +
  geom_boxplot()
