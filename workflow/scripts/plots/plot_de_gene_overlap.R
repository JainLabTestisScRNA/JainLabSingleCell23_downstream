library(tidyverse)
library(ggVennDiagram)
library(ComplexUpset)
de_objs <- read_rds("results/differential_expression/adult.df.germ_cell.reprocessed.de_results.rds")

de <- de_objs |>
  map_df(as_tibble,rownames="feature",.id="celltype") |>
  filter(FDR < 0.05) |>
  filter(str_detect(feature,"ENSMUSG")) |>
  group_by(celltype) |>
  summarise(de = list(feature)) |>
  deframe()


ggVennDiagram(de[c("3/Spermatocyte","4/Spermatocyte")])

ggVennDiagram(de[c("3/Spermatocyte","4/Spermatocyte","5/RoundSpermatid")])

ggVennDiagram(de[c("4/Spermatocyte","5/RoundSpermatid","6/RoundSpermatid","7/RoundSpermatid")])

ggVennDiagram(de[c("5/RoundSpermatid","6/RoundSpermatid","7/RoundSpermatid")])



