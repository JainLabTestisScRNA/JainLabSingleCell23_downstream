library(batchelor)
library(scran)
library(scuttle)
library(tidyverse)
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
# import data and gene variance calcs
# ------------------------------------------------------------------------------
#fls <- Sys.glob("results/single-cell-preproc/preprocessed/adult*.cellranger.sce.passing_cells.rds")
fls <- snakemake@input$sces

#dec_fls <- Sys.glob("results/single-cell-preproc/preprocessed/adult*.cellranger.dec.rds")
dec_fls <- snakemake@input$decs

names(fls) <- str_extract(fls, "adult\\d")
names(dec_fls) <- str_extract(dec_fls, "adult\\d")

sces <- map(fls, read_rds)
decs <- map(dec_fls, read_rds)[names(sces)]

# ------------------------------------------------------------------------------
# check if we need batch correction
# ------------------------------------------------------------------------------
# see https://bioconductor.org/books/3.17/OSCA.multisample/integrating-datasets.html#quick-start for motivations and general approach
# following the less convenient stepwise approach from here: https://bioconductor.org/books/3.17/OSCA.multisample/integrating-datasets.html#slower-setup
# in actuality. THis is to avoid redoing steps later that are very slow.
# will likely rename sample later, but want to retain batch info
sces <- map(sces, ~{.x$batch <- .x$Sample; .x})

# prep for combining
sces <- map(sces, ~{rowData(.x)$subset <- NULL;.x})
sces <- map(sces, ~{reducedDims(.x) <- NULL; .x})

# combine the sces
uncorrected <- combineCols(sces$adult1, sces$adult2, sces$adult3, sces$adult4)


# The normalized counts previously in here wouldn't be useful for between sample comparisons
# because they don't account for between batch differences in depth
uncorrected <- multiBatchNorm(uncorrected, batch = uncorrected$batch)

# still want hvgs made separately
combined.dec <- combineVar(decs)
hvg.var <- getTopHVGs(combined.dec, n=2500)

tes <- rowData(uncorrected)[rowData(uncorrected)$gene_biotype %in% "repeat_element",] |> rownames()
hvtes <- hvg.var[hvg.var %in% tes]
chosen <- hvg.var[!hvg.var %in% tes]
rowSubset(uncorrected) <- chosen
metadata(uncorrected)$highly.variable.tes <- hvtes
metadata(uncorrected)$highly.variable.genes <- chosen

# PCA - without influence of mito genes or TEs
uncorrected <- runPCA(uncorrected, subset_row=chosen)

g_pc_elbow_uncorrected <- tibble(percent.var = attr(reducedDim(uncorrected), "percentVar")) |>
  mutate(PC = row_number()) |>
  ggplot(aes(PC, percent.var)) +
  geom_point() +
  geom_vline(xintercept = ncol(reducedDim(uncorrected,"PCA")), color="red")

g_pca_uncorrected <- plotPCA(uncorrected, colour_by="batch")

# ------------------------------------------------------------------------------
# UMAP + TSNE uncorrected - spoiler - we def need batch correction
# ------------------------------------------------------------------------------

set.seed(2)
uncorrected <- runUMAP(uncorrected, dimred="PCA",pca=10)
set.seed(2)
uncorrected <- runTSNE(uncorrected, dimred="PCA",pca=10)

g_umap_uncorrected <- plotReducedDim(uncorrected, dimred = "UMAP", colour_by = "batch")
g_tsne_uncorrected <- plotReducedDim(uncorrected, dimred = "TSNE", colour_by = "batch")

# ------------------------------------------------------------------------------
# integration
# ------------------------------------------------------------------------------

# nuke the dimred stuff we added before
reducedDims(uncorrected) <- NULL
library(BiocParallel)

set.seed(1000101001)
mnn.out <- fastMNN(uncorrected, 
                   batch = uncorrected$batch, 
                   subset.row=chosen, 
                   correct.all = T, 
                   BSPARAM=BiocSingular::RandomParam(deferred = T),
                   BPPARAM = BiocParallel::MulticoreParam(8))

g_pca_corrected <- plotReducedDim(mnn.out, dimred = "corrected",colour_by = "batch")

# ------------------------------------------------------------------------------
# UMAP + TSNE corrected
# ------------------------------------------------------------------------------

set.seed(2)
mnn.out <- runUMAP(mnn.out, dimred="corrected")
set.seed(2)
mnn.out <- runTSNE(mnn.out, dimred="corrected")

g_umap_corrected <- plotReducedDim(mnn.out, dimred = "UMAP", colour_by = "batch")
g_tsne_corrected <- plotReducedDim(mnn.out, dimred = "TSNE", colour_by = "batch")


# ------------------------------------------------------------------------------
# complete the merging - get corrected/reconstructed values
# in the same sce as the original counts and metadata
# ------------------------------------------------------------------------------

stopifnot(all(rownames(mnn.out) == rownames(uncorrected)))
stopifnot(all(colnames(mnn.out) == colnames(uncorrected)))

corrected <- uncorrected

assay(corrected,"reconstructed") <- assay(mnn.out,"reconstructed")

reducedDims(corrected) <- reducedDims(mnn.out)

rowData(corrected)$rotation <- rowData(mnn.out)$rotation

metadata(corrected) <- metadata(mnn.out)

write_rds(corrected,snakemake@output$rds)
write_rds(g_pca_uncorrected, snakemake@output$g_pca_uncorrected)
write_rds(g_tsne_uncorrected, snakemake@output$g_tsne_uncorrected)
write_rds(g_umap_uncorrected, snakemake@output$g_umap_uncorrected)