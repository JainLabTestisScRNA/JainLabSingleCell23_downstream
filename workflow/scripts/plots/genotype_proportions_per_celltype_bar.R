Sys.setenv(R_PROFILE=".Rprofile")
source(Sys.getenv("R_PROFILE"))

library(tidyverse)
library(scater)
library(scran)
library(scuttle)
library(paletteer)

theme_set(theme_classic())

fl <- ifelse(exists("snakemake"),snakemake@input$sce,"results/integration/adult.sce.integrated.clustered.celltypes.rds")
sce <- read_rds(fl)

x <- makePerCellDF(sce) |>
  as_tibble(rownames="cell") |>
  dplyr::select(cell,label=celltype,genotype,batch)


g <- ggplot(x,aes(genotype,fill=label)) +
  geom_bar(position="fill") +
  theme(axis.text.x=element_text(angle=45,hjust=1)) +
  scale_fill_paletteer_d("RColorBrewer::Set3")

ggsave(snakemake@output$pdf,g)

g$data |> write_tsv(snakemake@output$tsv)