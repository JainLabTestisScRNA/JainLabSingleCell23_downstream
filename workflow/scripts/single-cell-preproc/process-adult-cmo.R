library(tidyverse)
library(DropletUtils)
library(scater)
library(scran)
library(scuttle)
library(ggdensity)
library(scDblFinder)
library(flexmix)
library(celda)
library(bluster)
library(PCAtools)

# ------------------------------------------------------------------------------
# files

sce_fl <- "results/cellranger-import/adult_9646_1.sce.assigned.rds"
sce_fl <- snakemake@input$sce

cmo_fl <- "data/cellranger/adult_9646_1/assignment_confidence_table.csv"
cmo_fl <- snakemake@input$sample_assignments

# ------------------------------------------------------------------------------
# import sce objs + add known cell types

sce <- read_rds(sce_fl)

colData(sce) <- as_tibble(colData(sce),rownames="cell") |>
  column_to_rownames("cell") |>
  DataFrame()

# ------------------------------------------------------------------------------
# doublet removal - mostly for any heterotypics unable to be removed
# via cellranger multi

# get logcounts
sce <- logNormCounts(sce,assay.type="counts")

# find doublets
set.seed(2)
scdbf <- scDblFinder(sce)

# get doublet finder data as tbl
dat_doublet <- colData(scdbf) %>% as_tibble(rownames = "cell")

metadata(sce)$doublets_removed <- dat_doublet |>
  filter(scDblFinder.class == "doublet") |>
  nrow()

# perform filt
sce <- sce[,scdbf$scDblFinder.class == "singlet"]

# ------------------------------------------------------------------------------
# basic dead/dying qc
is.mito <- grep("chrM",as.character(seqnames(sce)))

is.te <- rowData(sce)$gene_biotype %in% c("repeat_element") & !rowData(sce)$gene_biotype %in% c("7SLRNA")

sce <- addPerCellQCMetrics(sce,subsets=(list(MT=is.mito,
                                             TE=is.te,
                                             `7SLRNA`=names(sce)[str_detect(names(sce),"7SL")])))

discard.mito.madTop <- isOutlier(sce$subsets_MT_percent, type="higher", nmads=3)
discard.mito.mad1 <- isOutlier(sce$subsets_MT_percent, type="higher", nmads=2)

g_mito_cuts <- tibble(total.decont.umis = sce$sum,
       mito.percent = sce$subsets_MT_percent) |>
  ggplot(aes(total.decont.umis,mito.percent)) +
  geom_hline(yintercept = attr(discard.mito.madTop,"thresholds")["higher"], color="red") +
  geom_hline(yintercept = attr(discard.mito.mad1,"thresholds")["higher"], color="blue") +
  geom_point() +
  scale_x_log10()

sce$mito_warning <- as.vector(discard.mito.mad1)
sce <- sce[,!discard.mito.madTop]

metadata(sce)$n_mito_discarded <- sum(discard.mito.madTop)
metadata(sce)$n_mito_warning <- sum(sce$mito_warning)

# ------------------------------------------------------------------------------
# remove low dimension cells - very few features expressed, likely another class of dying cells
# but either way not super informative and tends to screw up trajectory analyses later on

sce$nexprs <- nexprs(sce) 

discard.nexpres.madlower <- isOutlier(sce$nexprs, type="lower", nmads=3,log = T)

sce <- sce[,!discard.nexpres.madlower]

metadata(sce)$n_too_few_genes_exprs <- sum(discard.nexpres.madlower)


# ------------------------------------------------------------------------------
# renormalize (previously used simple norm for scbdlfinder), deconvolution norm for actual analysis

sce$sizeFactor <- NULL
logcounts(sce) <- NULL

clust.pre <- quickCluster(sce)

# occasionally get warning about negative size factors - use positive=T to enforce positivity
sce <- computeSumFactors(sce,cluster=clust.pre, min.mean=0.2, positive=T)

sce <- logNormCounts(sce)

# ------------------------------------------------------------------------------
# variable feature selection

dec <- modelGeneVar(sce)
hvg.var <- getTopHVGs(dec, prop = 0.1)

tes <- rowData(sce)[rowData(sce)$gene_biotype %in% "repeat_element",] |> rownames()

hvtes <- hvg.var[hvg.var %in% tes]

#& 
chosen <- hvg.var[!hvg.var %in% tes & !hvg.var %in% subset(rowData(sce),grepl("^mt-",gene_name))$gene_name]

rowSubset(sce) <- chosen
metadata(sce)$highly.variable.tes <- hvtes
metadata(sce)$highly.variable.genes <- chosen

# ------------------------------------------------------------------------------
# PCA - without influence of mito genes or TEs
sce <- fixedPCA(sce, subset.row=chosen)

percent.var <- attr(reducedDim(sce,"PCA"), "percentVar")

chosen.elbow <- findElbowPoint(percent.var)

g_pc_elbow <- tibble(percent.var = attr(reducedDim(sce), "percentVar")) |>
  mutate(PC = row_number()) |>
  ggplot(aes(PC, percent.var)) +
  geom_point() +
  geom_vline(xintercept = chosen.elbow, color="red")

# ------------------------------------------------------------------------------
# UMAP
set.seed(2)
sce <- runUMAP(sce, dimred="PCA")
set.seed(2)
sce <- runTSNE(sce, dimred="PCA")

# ------------------------------------------------------------------------------
# export
write_rds(sce, snakemake@output$rds)
write_rds(g_mito_cuts, snakemake@output$g_mito_cuts)
write_rds(g_pc_elbow, snakemake@output$g_pc_elbow)