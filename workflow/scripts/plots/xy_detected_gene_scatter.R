Sys.setenv(R_PROFILE=".Rprofile")
source(Sys.getenv("R_PROFILE"))

library(tidyverse)
library(scater)
library(scran)
library(scuttle)

sce <- read_rds(ifelse(exists("snakemake"),snakemake@input$sce,
                       "results/germ_cells/adult.sce.germ_cell.both_genotypes.subclustered.reintegrated.rds"))

theme_set(theme_classic())

dat <- colData(sce) |>
  as_tibble(rownames="cell") |>
  dplyr::select(cell,genotype,label,subsets_TE_percent,subsets_TE_sum,total,nexprs,matches("chr")) |>
  filter(str_detect(label,"Spermatid|Spermatocyte")) |>
  pivot_longer(-c(cell,genotype,label,contains("subsets_TE"),total,nexprs,subsets_chrX_detected))

plot_vs_chrX <- function(chr="chrY") {
  message(chr)
  dat |>
    filter(name==sprintf("subsets_%s_detected",chr)) |>
    ggplot(aes(subsets_chrX_detected,value)) +
    geom_point(size=rel(0.5)) +
    facet_wrap(genotype~label,nrow=2,scales="free") +
    xlab("chrX genes detected") +
    ylab(sprintf("%s genes detected",chr)) +
    theme(aspect.ratio = 1)
}

pdf(snakemake@output$pdf, width = 8.5, height = 6)
rowRanges(sce) |>
  as_tibble() |>
  dplyr::select(seqnames) |>
  distinct() |>
  filter(sprintf("subsets_%s_detected",seqnames) %in% dat$name) |>
  pull(seqnames) |>
  droplevels() |>
  as.character() |>
  sort(decreasing = T) |>
  map(plot_vs_chrX)
dev.off()


write_tsv(dat,snakemake@output$tsv)