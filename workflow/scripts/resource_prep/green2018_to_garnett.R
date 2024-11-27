Sys.setenv(R_PROFILE=".Rprofile")
source(Sys.getenv("R_PROFILE"))

library(tidyverse)
library(readxl)

# https://doi.org/10.1016/j.devcel.2018.07.025
x <- read_xlsx("data/green2018_supp2.xlsx",range = "A4:H4912")


x <- x |> filter(CellType != "Unknown")

x |>
  mutate(p_val = as.numeric(p_val)) |>
  filter(p_val <= 0.0001) |>
  group_by(CellType) |> 
  slice_max(avg_diff, n=100,with_ties = T) |>
  summarize(expressed = paste(gene,collapse=", ")) |>
  mutate(ix= row_number()) |>
  mutate(CellType=paste0(">",CellType)) |>
  mutate(expressed=paste0("expressed: ",expressed)) |>
  mutate(spacer = "\n") |>
  pivot_longer(c(CellType,expressed,spacer),names_to = "entry") |>
  dplyr::select(value) |>
  write_tsv(file = "data/green2018_garnett_adult_testis_classifier.txt", col_names = F)

