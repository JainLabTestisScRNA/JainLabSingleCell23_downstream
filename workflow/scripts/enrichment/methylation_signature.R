Sys.setenv(R_PROFILE=".Rprofile")
source(Sys.getenv("R_PROFILE"))

library(tidyverse)
library(scater)
library(scran)
library(scuttle)
library(fgsea)
library(GSEABase)
library(AUCell)
library(doMC)

# ------------------------------------------------------------------------------
# pick signatures
# ------------------------------------------------------------------------------
gs <- read_rds("results/enrichment/gsea_prerank.rds")

gs <- map_df(gs,as_tibble,.id="label") |>
  filter(str_detect(label,"Spermatogonia|Spermatocyte")) |>
  filter(str_detect(pathway,"GOBP")) |>
  group_by(label, sign(NES)) |>
  slice_min(pval,n=10) |>
  pull(pathway) |>
  unique()

gmt0 <- GSEABase::getGmt(ifelse(exists("snakemake"),snakemake@input$msigdb,"data/m5.all.v2023.2.Mm.symbols.gmt"))

gmt <- gmt0[str_detect(names(gmt0),"GOCC")]
gmt <- gmt[map_lgl(gmt,~{dplyr::between(length(.x@geneIds),25,100)})]

sce <- read_rds(ifelse(exists("snakemake"),snakemake@input$sce,"results/germ_cells/adult.sce.germ_cell.both_genotypes.subclustered.reintegrated.rds"))

mat <- counts(sce)
rownames(mat) <- rowData(sce)$gene_name

rankings <- AUCell_buildRankings(mat,splitByBlocks = T,
                                 plotStats=F, verbose=FALSE)

# default aucMaxRank = 0.05
scsig.aucs <- AUCell_calcAUC(gmt, rankings,aucMaxRank = ceiling(0.1 * nrow(rankings)),verbose = T,nCores = 4)

scsig.aucs <- t(assay(scsig.aucs))

stopifnot(all(rownames(scsig.aucs) == colnames(sce)))


dat <- scsig.aucs |>
  as_tibble(rownames="cell") |>
  mutate(label = sce$label, genotype = sce$genotype) |>
  dplyr::relocate(label,genotype,.after=cell) |>
  pivot_longer(-c(cell,label,genotype))


dat2 <- dat |>
  filter(str_detect(label,"Spermatogonia|Spermatocyte")) |>
  group_by(label,genotype,name) |>
  summarise(scores = list(value)) |>
  pivot_wider(names_from = genotype,values_from = scores) |>
  mutate(lfc = map2_dbl(MUT,WT,~log2((median(.x)+0.001)/(median(.y)+0.001)))) |>
  mutate(wilc = map2(MUT,WT,~wilcox.test(.x,.y))) |>
  mutate(tid = map(wilc,broom::tidy)) |>
  unnest(tid)

dat2 |>
  #arrange(-p.value) |>
  mutate(padj = p.adjust(p.value,method="BH")) |>
  filter(padj < 0.001) |>
  arrange(lfc) |>
  pull(name) |>
  tail(20) |>
  map(plot_signature)













plot_signature <- function(x) {
  return(dat |>
           filter(name==x) |>
           ggplot(aes(label,value,fill=genotype)) +
           facet_wrap(~name, scales="free") +
           geom_boxplot())
}

map(names(gmt),plot_signature)
