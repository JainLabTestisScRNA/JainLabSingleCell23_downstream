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

# ------------------------------------------------------------------------------
# files

#sce_fl <- "results/single-cell-preproc/count/import/adult2.cellranger.sce.rds"
sce_fl <- snakemake@input$sce

#raw_fl <- "results/single-cell-preproc/count/raw/adult2.cellranger.sce.rds"
raw_fl <- snakemake@input$raw
# ------------------------------------------------------------------------------

# import sce objs 

sce <- read_rds(sce_fl)
raw <- read_rds(raw_fl)

colData(sce) <- as_tibble(colData(sce),rownames="cell") |>
  column_to_rownames("cell") |>
  DataFrame()


#-------------------------------------------------------------------------------
# ambient rna removal

# set.seed(2)
# dcx.sce <- decontX(sce, seed=12345, background=raw[,!raw$Barcode %in% sce$Barcode])
# 
# counts(sce) <- round(decontXcounts(dcx.sce))
# 
# g_dcx <- celda::plotDecontXContamination(dcx.sce)

# ------------------------------------------------------------------------------
# filtering after ambient removal

bcrank <- barcodeRanks(counts(sce))

g_kneeplot <- as_tibble(bcrank) |>
  ggplot(aes(rank,total)) +
  geom_line() +
  scale_y_log10() +
  scale_x_log10() +
  geom_hline(yintercept = metadata(bcrank)$inflection, color="red") +
  geom_hline(yintercept = metadata(bcrank)$knee, color="blue") +
  geom_hline(yintercept = 1000, color="darkgreen")


sce <- sce[,colSums(counts(sce)) > 500]

# ------------------------------------------------------------------------------
# doublet removal

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

discard.mito.madTop <- isOutlier(sce$subsets_MT_percent, type="higher", nmads=5)
mito_cut <- attr(discard.mito.madTop,"thresholds")["higher"]
#mito_cut <- 10
#discard.mito.madTop <- sce$subsets_MT_percent > mito_cut

g_mito_cuts <- tibble(total.decont.umis = sce$sum,
       mito.percent = sce$subsets_MT_percent) |>
  ggplot(aes(total.decont.umis,mito.percent)) +
  geom_hline(yintercept = mito_cut, color="red") +
  geom_point() +
  scale_x_log10()

sce <- sce[,!discard.mito.madTop]

metadata(sce)$n_mito_discarded <- sum(discard.mito.madTop)

# ------------------------------------------------------------------------------
# remove low dimension cells - very few features expressed, likely another class of dying cells
# but either way not super informative and tends to screw up trajectory analyses later on

sce$nexprs <- nexprs(sce) 

discard.nexpres.madlower <- isOutlier(sce$nexprs, type="lower", nmads=3)

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
dec.tech <- modelGeneVarByPoisson(sce)
sce <- denoisePCA(sce, dec.tech, subset.row=chosen)

g_pc_elbow <- tibble(percent.var = attr(reducedDim(sce), "percentVar")) |>
  mutate(PC = row_number()) |>
  ggplot(aes(PC, percent.var)) +
  geom_point() +
  geom_vline(xintercept = ncol(reducedDim(sce,"PCA")), color="red")


# ------------------------------------------------------------------------------
# UMAP

set.seed(2)
sce <- runUMAP(sce, dimred="PCA")
set.seed(2)
sce <- runTSNE(sce, dimred="PCA")

#plotReducedDim(sce[,order(sce$subsets_MT_percent)], dimred = "UMAP", colour_by = "subsets_MT_percent")
#plotReducedDim(sce[,order(sce$subsets_TE_percent)], dimred = "UMAP", colour_by = "subsets_TE_percent")
#plotReducedDim(sce, dimred = "UMAP", colour_by = "subsets_7SLRNA_percent")
#plotReducedDim(sce[,order(sce$nexprs)], dimred = "UMAP", colour_by = "nexprs")
#plotReducedDim(sce[,order(sce$detected)], dimred = "UMAP", colour_by = "detected")
#plotReducedDim(sce[,order(sce$sum)], dimred = "UMAP", colour_by = "sum")
#plotReducedDim(sce[,order(sce$total)], dimred = "UMAP", colour_by = "total")

#tibble(x = sce$sum) |>
#  ggplot(aes(x)) +
#  geom_density() +
#  scale_x_log10()
# ------------------------------------------------------------------------------
# skipping clustering here, bc this will be integrated and reclustered later

# ------------------------------------------------------------------------------
# export
write_rds(sce, snakemake@output$rds)
write_rds(dec, snakemake@output$dec)
write_rds(NULL, snakemake@output$g_dcx)
write_rds(g_mito_cuts, snakemake@output$g_mito_cuts)
write_rds(g_pc_elbow, snakemake@output$g_pc_elbow)
write_rds(g_kneeplot, snakemake@output$g_kneeplot)
