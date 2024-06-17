Sys.setenv(R_PROFILE=".Rprofile")
source(Sys.getenv("R_PROFILE"))

library(tidyverse)
library(corrr)
dat <- read_tsv(ifelse(exists("snakemake"),
                       snakemake@input$tsv,
                       "results/differential_expression/adult.tbl.germ_cell.reprocessed.de.tsv.gz"))

mat <- dat |>
  dplyr::select(celltype,feature,logFC) |>
  pivot_wider(names_from = celltype,values_from = logFC)

nms <- colnames(mat)[2:ncol(mat)]

g <- mat |> correlate() |> corrr::stretch() |>
  drop_na() |>
  ggplot(aes(x,y,color=r)) +
  geom_point(size=rel(5)) +
  scale_color_distiller(palette = 2,direction = 1) +
  theme_classic() +
  theme(axis.text.x = element_text(angle=45,hjust=1)) +
  xlab("") +
  ylab("")


pdf(snakemake@output$pdf)

g
dev.off()

write_tsv(g$data,snakemake@output$tsv)
