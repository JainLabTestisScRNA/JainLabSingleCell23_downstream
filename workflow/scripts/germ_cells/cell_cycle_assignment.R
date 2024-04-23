library(tidyverse)
library(scater)
library(scran)
library(scuttle)
library(AUCell)
library(GSEABase)

gmt <- getGmt("data/m5.all.v2023.2.Mm.symbols.gmt")

meio <- gmt[str_detect(names(gmt),"MEIO")]

names(meio)

#sigs <- read_tsv("results/resources/cell_cycle_genes.tsv")
#sigs <- sigs |>
#  dplyr::select(name,ENSEMBL) |>
#  drop_na() |>
#  distinct()
#sigs <- split(sigs$ENSEMBL,sigs$name)


sce <- read_rds("results/germ_cells/adult.sce.germ_cell.both_genotypes.subclustered.reintegrated.rds")

sce <- sce[,sce$genotype == "WT"]

mat <- counts(sce)
rownames(mat) <- rowData(sce)$gene_name


rankings <- AUCell_buildRankings(mat,splitByBlocks = T,
                                        plotStats=T, verbose=FALSE)

# default aucMaxRank = 0.05
scsig.aucs <- AUCell_calcAUC(meio, rankings,aucMaxRank = ceiling(0.05 * nrow(rankings)))

scsig.aucs <- t(assay(scsig.aucs))

stopifnot(all(rownames(scsig.aucs) == colnames(sce)))

names(meio)

dat <- scsig.aucs |>
  as_tibble(rownames="cell") |>
  mutate(label = sce$label) |>
  dplyr::relocate(label,.after=cell)

plot_meio_signature <- function(x) {
  return(plotColData(sce,x="label",y=I(scsig.aucs[,x]),other_fields = "genotype") + 
    facet_wrap(~genotype) + 
    theme(axis.text.x = element_text(angle=45, hjust=1)) +
    labs(title=x))
}

pdf(snakemake@output$pdf,height = 3, width = 6)
names(meio) |>
  set_names() |>
  map(plot_meio_signature)
dev.off()

write_tsv(dat,snakemake@output$tsv)

