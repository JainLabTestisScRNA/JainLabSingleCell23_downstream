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
library(readxl)

# ------------------------------------------------------------------------------
# grr signature (hill et al 2018)
# ------------------------------------------------------------------------------
grr <- read_xlsx("~/Downloads/41586_2018_BFnature25964_MOESM3_ESM.xlsx",col_names = "gene_name") |> pull(gene_name)

sce <- read_rds(ifelse(exists("snakemake"),snakemake@input$sce,"results/germ_cells/adult.sce.germ_cell.both_genotypes.subclustered.reintegrated.rds"))

mat <- counts(sce)
rownames(mat) <- rowData(sce)$gene_name

rankings <- AUCell_buildRankings(mat,splitByBlocks = T,
                                 plotStats=F, verbose=FALSE)

# default aucMaxRank = 0.05, but this is a very small list ~45 genes, so doesn't make sense to
scsig.aucs <- AUCell_calcAUC(grr, rankings,aucMaxRank = ceiling(1 * nrow(rankings)),verbose = T,nCores = 4)

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


plot_signature <- function(x) {
  return(dat |>
           filter(name==x) |>
           ggplot(aes(label,value,fill=genotype)) +
           facet_wrap(~name, scales="free") +
           geom_boxplot())
}

dat2 |>
  #arrange(-p.value) |>
  mutate(padj = p.adjust(p.value,method="BH")) |>
  filter(padj < 0.05) |>
  arrange(lfc) |>
  pull(name) |>
  tail(20) |>
  map(plot_signature)


plotReducedDim(sce,dimred="corrected",ncomponents = c(1,3),colour_by = I(scsig.aucs[,1]),other_fields = "genotype") + facet_wrap("genotype")
plotReducedDim(sce,dimred="corrected",ncomponents = c(1,3),colour_by = "label",other_fields = "genotype") + facet_wrap("genotype")


plotExpression(sce,"Pou5f1",x="label",swap_rownames = "gene_name",other_fields = "genotype") + facet_wrap("genotype")

plotExpression(sce,"Pou5f1",x="L1MdA_I_3end",swap_rownames = "gene_name",colour_by = "label",other_fields = c("genotype","label")) +
  facet_wrap(~genotype + label)

# ------------------------------------------------------------------------------
# grr in de
# ------------------------------------------------------------------------------
de <- read_tsv("results/differential_expression/adult.tbl.germ_cell.reprocessed.de.tsv.gz")

lkup <- dplyr::select(as_tibble(rowData(sce)), feature=gene_id,gene_name)

de <- left_join(de,lkup)


de |> filter(gene_name %in% grr) |> filter(FDR < 0.1) |> print(n=Inf)

# ------------------------------------------------------------------------------
# grr in de
# ------------------------------------------------------------------------------