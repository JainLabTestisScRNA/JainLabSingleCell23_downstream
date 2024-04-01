library(tidyverse)
library(readxl)

# https://doi.org/10.1016/j.devcel.2018.07.025
x <- read_xlsx("data/green2018_supp3.xlsx",range = "A6:H2625",sheet = "A.GC_Markers")

x |>
  mutate(p_val = as.numeric(p_val)) |>
  filter(p_val <= 0.01) |>
  mutate(GermCellState = paste0("GC",GermCellState)) |>
  group_by(GermCellState) |>
  summarize(expressed = paste(gene,collapse=", ")) |>
  mutate(ix= row_number()) |>
  mutate(GermCellState=paste0(">",GermCellState)) |>
  mutate(expressed=paste0("expressed: ",expressed)) |>
  mutate(spacer = "\n") |>
  pivot_longer(c(GermCellState,expressed,spacer),names_to = "entry") |>
  dplyr::select(value) |>
  write_tsv(file = "data/green2018_garnett_adult_germcell_classifier.txt", col_names = F)

