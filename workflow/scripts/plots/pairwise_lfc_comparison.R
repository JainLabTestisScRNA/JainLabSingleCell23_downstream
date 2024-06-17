Sys.setenv(R_PROFILE=".Rprofile")
source(Sys.getenv("R_PROFILE"))

library(tidyverse)

dat <- read_tsv(ifelse(exists("snakemake"),
                              snakemake@input$tsv,
                              "results/differential_expression/adult.tbl.germ_cell.reprocessed.de.tsv.gz"))

# slashes cause problems with aes_string
#dat$celltype <- dat$celltype |> str_replace_all("/",".")

mat <- dat |>
  dplyr::select(celltype,feature,logFC) |>
  pivot_wider(names_from = celltype,values_from = logFC)

nms <- colnames(mat)[2:ncol(mat)]

plot_df <- combn(nms,2,simplify = F) |>
  enframe() |>
  unnest(value) |>
  group_by(name) |>
  mutate(axis=c("x","y")[row_number()]) |>
  pivot_wider(names_from = axis,values_from = value) |>
  ungroup()



z <- function(X,Y) {
  ggplot(mat,aes(.data[[X]],.data[[Y]])) +
    ggdensity::geom_hdr_points() +
    #geom_point() +
    ggpubr::stat_cor() +
    theme(aspect.ratio=1)
}

plot_df <- plot_df |>
  mutate(gg=map2(x,y,z))

pdf(snakemake@output$pdf)
plot_df$gg
dev.off()

write_tsv(mat,snakemake@output$tsv)