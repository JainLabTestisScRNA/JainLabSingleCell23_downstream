Sys.setenv(R_PROFILE=".Rprofile")
source(Sys.getenv("R_PROFILE"))

library(tidyverse)
library(readxl)

# https://doi.org/10.1016/j.devcel.2018.07.025
x <- read_xlsx("data/green2018_supp3.xlsx",range = "A6:H2625",sheet = "A.GC_Markers")

x |>
  filter(!str_detect(gene,"mt-|Rik|^RP|^Rps")) |>
  mutate(p_val = as.numeric(p_val)) |>
  mutate(GermCellState = paste0("GC",GermCellState)) |>
  group_by(GermCellState) |>
  slice_max(pct.1-pct.2,n=20) |>
  summarize(expressed = paste(gene,collapse=", ")) |>
  mutate(ix= row_number()) |>
  mutate(GermCellState=paste0(">",GermCellState)) |>
  mutate(expressed=paste0("expressed: ",expressed)) |>
  mutate(spacer = "\n") |>
  pivot_longer(c(GermCellState,expressed,spacer),names_to = "entry") |>
  dplyr::select(value) |>
  write_tsv(file = "data/green2018_garnett_adult_germcell_classifier.txt", col_names = F)

