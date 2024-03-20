library(batchelor)
library(scran)
library(scuttle)
library(tidyverse)
library(DropletUtils)
library(scater)
library(ggdensity)
library(scDblFinder)
library(flexmix)
library(celda)
library(bluster)
library(monocle3)
library(garnett)
library(org.Mm.eg.db)


# ------------------------------------------------------------------------------
# import data 
# ------------------------------------------------------------------------------
fls <- Sys.glob("results/single-cell-preproc/preprocessed/adult_96*.cellranger.sce.passing_cells.rds")
fls <- snakemake@input$sces
names(fls) <- str_extract(fls, "adult_\\d+_\\d")

sces <- map(fls, read_rds)

bc_lists <- sces |> 
  map_df(~{
    makePerCellDF(.x) |>
      as_tibble() |>
      dplyr::select(Barcode)
  },.id="run") |>
  group_by(run) |>
  summarise(Barcode = list(Barcode)) |>
  deframe()

library(ggVennDiagram)
gg_barcode_overlap <- ggVennDiagram(bc_lists)


# rename to be safe
colnames(sces$adult_9646_2) <- colnames(sces$adult_9646_2) |> str_replace("-1","-2")
colnames(sces$adult_9682_1) <- colnames(sces$adult_9682_1) |> str_replace("-1","-3")
colnames(sces$adult_9682_2) <- colnames(sces$adult_9682_2) |> str_replace("-1","-4")

# per email from Devanshi 
# 9646 wt is cmo 307, 9646 mut is cmo 308
# 9682 wt is cmo 305, 9682 mut is cmo 306
sces$adult_9646_1$genotype <- sces$adult_9646_1$Assignment |> magrittr::equals("CMO307") |> if_else("WT","MUT")
sces$adult_9646_2$genotype <- sces$adult_9646_2$Assignment |> magrittr::equals("CMO307") |> if_else("WT","MUT")
sces$adult_9682_1$genotype <- sces$adult_9682_1$Assignment |> magrittr::equals("CMO305") |> if_else("WT","MUT")
sces$adult_9682_2$genotype <- sces$adult_9682_2$Assignment |> magrittr::equals("CMO305") |> if_else("WT","MUT")


# ------------------------------------------------------------------------------
# set up for merging
# ------------------------------------------------------------------------------
# see https://bioconductor.org/books/3.17/OSCA.multisample/integrating-datasets.html#quick-start for motivations and general approach
# following the less convenient stepwise approach from here: https://bioconductor.org/books/3.17/OSCA.multisample/integrating-datasets.html#slower-setup
# in actuality. THis is to avoid redoing steps later that are very slow.
# will likely rename sample later, but want to retain batch info

# assign batches based on devanshi's 2/6/24 email titled `Sc RNA seq genotypes`
sces <- imap(sces, ~{.x$batch <- paste0("batch",str_extract(.y,"\\d{4}_\\d")); .x})

# make sure cell names are unique - original barcodes - <10 shared between batches - are preserved in coldata
#colnames(sces$adult_9646_combined) <- colnames(sces$adult_9646_combined) |> str_replace("-1","-2")

# prep for combining
sces <- map(sces, ~{rowData(.x)$subset <- NULL;.x})
sces <- map(sces, ~{reducedDims(.x) <- NULL; .x})
sces <- map(sces, ~{sizeFactors(.x) <- NULL; .x})

# these cols not shared between batches
sces$adult_9646_1$CMO307 <- NULL
sces$adult_9646_1$CMO308 <- NULL
sces$adult_9646_2$CMO307 <- NULL
sces$adult_9646_2$CMO308 <- NULL

sces$adult_9682_1$CMO305 <- NULL
sces$adult_9682_1$CMO306 <- NULL
sces$adult_9682_2$CMO305 <- NULL
sces$adult_9682_2$CMO306 <- NULL

# make sure all obj have the same genes.
# because these were quantified via the same pipeline
# usually unnecessary unless something weird happened
# or we pre-filtered genes
universe <- map(sces,rownames) |> Reduce(intersect,x=_)
sces <- map(sces,~{.x[universe,]})

# ------------------------------------------------------------------------------
# check if batch correction needed
# ------------------------------------------------------------------------------
# combine the sces
uncorrected0 <- cbind(sces$adult_9646_1,sces$adult_9646_2,sces$adult_9682_1,sces$adult_9682_2)
metadata(uncorrected0) <- list()

# The normalized counts previously in here wouldn't be useful for between sample comparisons
# because they don't account for between batch differences in depth
uncorrected0 <- multiBatchNorm(uncorrected0, batch = uncorrected0$batch)

# get hvgs
decs <- map(sces, modelGeneVar)
combined.dec <- combineVar(decs)
hvg.var <- getTopHVGs(combined.dec, n=10000)

tes <- rowData(uncorrected0)[rowData(uncorrected0)$gene_biotype %in% "repeat_element",] |> rownames()
hvtes <- hvg.var[hvg.var %in% tes]
chosen <- hvg.var[!hvg.var %in% tes]
rowSubset(uncorrected0) <- chosen
metadata(uncorrected0)$highly.variable.tes <- hvtes
metadata(uncorrected0)$highly.variable.genes <- chosen

# ------------------------------------------------------------------------------
# UMAP + TSNE uncorrected - we may not actually need batch correction with current
# preprocesing approach - only TSNE vis seems to have clusters that don't overlap
# ------------------------------------------------------------------------------
library(PCAtools)
# PCA - without influence of mito genes or TEs
uncorrected <- runPCA(uncorrected0, subset_row=chosen, ncomponents=50)

percent.var <- attr(reducedDim(uncorrected,"PCA"), "percentVar")

chosen.elbow <- findElbowPoint(percent.var)

g_pc_elbow_uncorrected <- tibble(percent.var = attr(reducedDim(uncorrected), "percentVar")) |>
  mutate(PC = row_number()) |>
  ggplot(aes(PC, percent.var)) +
  geom_point() +
  geom_vline(xintercept = chosen.elbow, color="red")

metadata(uncorrected)$pca_elbow <- chosen.elbow

set.seed(2)
uncorrected <- runUMAP(uncorrected, dimred="PCA",pca=chosen.elbow)


# ------------------------------------------------------------------------------
# integration
# https://bioconductor.org/books/3.17/OSCA.multisample/integrating-datasets.html#mnn-correction
# "...k... can be interpreted as the minimum anticipated frequency of any shared cell type or state in each batch."
# ------------------------------------------------------------------------------
set.seed(2)
mnn.out <- fastMNN(uncorrected0, k=10,d = 50,
                   batch = uncorrected0$batch, 
                   subset.row=chosen, # no tes involved
                   correct.all = T, 
                   BSPARAM=BiocSingular::RandomParam(deferred = T))

corrected <- uncorrected0
reducedDims(corrected) <- reducedDims(mnn.out)
assay(corrected,"reconstructed") <- assay(mnn.out,"reconstructed")

# ------------------------------------------------------------------------------
# UMAP
# ------------------------------------------------------------------------------

set.seed(2)
corrected <- runUMAP(corrected, dimred="corrected", pca=chosen.elbow)

write_rds(corrected,snakemake@output$rds)
write_rds(uncorrected, snakemake@output$uncorrected)
write_rds(mnn.out, snakemake@output$mnn_out)