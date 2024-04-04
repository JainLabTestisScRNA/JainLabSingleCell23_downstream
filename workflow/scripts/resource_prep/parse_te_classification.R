library(tidyverse)

raw <- read_delim("data/dfam_classification.txt",col_names = c("dfam_id","dfam_name","classification_long","len","classification"))

# extract name from single quotes (and colon)
df <- raw |>
  mutate(dfam_name = str_extract(dfam_name,"(?<=').+(?=':)"))

# extract len
df <- df |> mutate(len =as.integer(str_extract(len,"(?<==).+")))

df <- dplyr::relocate(df,classification,.after="dfam_name")

write_tsv(df,"data/dfam_classification.parsed.txt")
