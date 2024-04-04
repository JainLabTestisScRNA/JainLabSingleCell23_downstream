library(tidyverse)
library(scuttle)
library(scran)
library(scater)
library(patchwork)

theme_set(theme_classic())

fl <- ifelse(exists("snakemake"),
             snakemake@input$sce,
             "results/integration/adult.sce.integrated.clustered.celltypes.rds")

sce <- read_rds(fl)

# this script plots a diagnostic suggested by Aaron lun for
# assessing issues with deconvolution-based size factors
#https://support.bioconductor.org/p/95017/
# https://bioconductor.org/books/3.18/OSCA.basic/normalization.html

x <- colSums(counts(sce)) |>
  enframe(name = "cell",value = "library_size")

x$deconvSizeFactors <- sizeFactors(sce)

x$librarySizeFactors <- librarySizeFactors(sce)

g <- ggplot(x,aes(library_size,deconvSizeFactors)) +
  geom_point()

g2 <- ggplot(x,aes(library_size,librarySizeFactors)) +
  geom_point()

g3 <- ggplot(x,aes(librarySizeFactors,deconvSizeFactors)) +
  geom_point()

gx <- g + g2 + g3

write_tsv(x, snakemake@output$tsv)
ggsave(snakemake@output$pdf,gx)