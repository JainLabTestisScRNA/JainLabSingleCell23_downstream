Sys.setenv(R_PROFILE=".Rprofile")
source(Sys.getenv("R_PROFILE"))

library(tidyverse)
library(scater)
library(scran)
library(scuttle)
library(AUCell)
library(GSEABase)

gmt <- getGmt("data/m5.all.v2023.2.Mm.symbols.gmt")

sigs <- gmt[str_detect(names(gmt),"MEIO|PIRNA|TRANSPOSIT|SPERM|REPLICATION|DNA_METHY")]
sce <- read_rds("results/germ_cells/adult.sce.germ_cell.both_genotypes.subclustered.reintegrated.rds")

mat <- counts(sce)

rownames(mat) <- rowData(sce)$gene_name

rankings <- AUCell_buildRankings(mat,splitByBlocks = T,
                                 plotStats=T, verbose=FALSE)

# default aucMaxRank = 0.05
scsig.aucs <- AUCell_calcAUC(sigs, rankings,aucMaxRank = ceiling(0.1 * nrow(rankings)))

scsig.aucs <- t(assay(scsig.aucs))

stopifnot(all(rownames(scsig.aucs) == colnames(sce)))

dat <- scsig.aucs |>
  as_tibble(rownames="cell") |>
  mutate(label = sce$label,
         genotype = sce$genotype) |>
  dplyr::relocate(label,genotype,.after=cell)


dat

plot_signature <- function(x) {
  d <- dplyr::select(dat,cell,label,genotype,matches(x))
  
  ggplot(d,aes_string("label",x,fill="genotype")) +
    geom_boxplot() +
    theme(axis.text.x = element_text(angle=45,hjust=1))
}

pdf(snakemake@output$pdf,height = 5, width = 6)
names(sigs) |>
  set_names() |>
  map(plot_signature)
dev.off()

write_tsv(dat,snakemake@output$tsv)

