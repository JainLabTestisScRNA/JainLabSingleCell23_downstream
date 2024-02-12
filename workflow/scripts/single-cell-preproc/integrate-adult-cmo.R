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
fls <- Sys.glob("results/single-cell-preproc/preprocessed/adult_96*combined.cellranger.sce.passing_cells.rds")
fls <- snakemake@input$sces
names(fls) <- str_extract(fls, "adult_\\d+_combined")

sces <- map(fls, read_rds)

# per email from Devanshi 
# 9646 wt is cmo 307, 9646 mut is cmo 308
# 9682 wt is cmo 305, 9682 mut is cmo 306
sces$adult_9646_combined$genotype <- sces$adult_9646_combined$Assignment |> magrittr::equals("CMO307") |> if_else("WT","MUT")
sces$adult_9682_combined$genotype <- sces$adult_9682_combined$Assignment |> magrittr::equals("CMO305") |> if_else("WT","MUT")

# ------------------------------------------------------------------------------
# set up for merging
# ------------------------------------------------------------------------------
# see https://bioconductor.org/books/3.17/OSCA.multisample/integrating-datasets.html#quick-start for motivations and general approach
# following the less convenient stepwise approach from here: https://bioconductor.org/books/3.17/OSCA.multisample/integrating-datasets.html#slower-setup
# in actuality. THis is to avoid redoing steps later that are very slow.
# will likely rename sample later, but want to retain batch info

# assign batches based on devanshi's 2/6/24 email titled `Sc RNA seq genotypes`
sces <- imap(sces, ~{.x$batch <- paste0("batch",str_extract(.y,"(?<=_)\\d+(?=_combined)")); .x})

# make sure cell names are unique - original barcodes - <10 shared between batches - are preserved in coldata
colnames(sces$adult_9646_combined) <- colnames(sces$adult_9646_combined) |> str_replace("-1","-2")

# prep for combining
sces <- map(sces, ~{rowData(.x)$subset <- NULL;.x})
sces <- map(sces, ~{reducedDims(.x) <- NULL; .x})
sces <- map(sces, ~{sizeFactors(.x) <- NULL; .x})

# these cols not shared between batches
sces$adult_9646_combined$CMO307 <- NULL
sces$adult_9646_combined$CMO308 <- NULL
sces$adult_9682_combined$CMO305 <- NULL
sces$adult_9682_combined$CMO306 <- NULL

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
uncorrected0 <- cbind(sces$adult_9646_combined,sces$adult_9682_combined)
metadata(uncorrected0) <- list()

# The normalized counts previously in here wouldn't be useful for between sample comparisons
# because they don't account for between batch differences in depth
uncorrected0 <- multiBatchNorm(uncorrected0, batch = uncorrected0$batch)

# get hvgs
decs <- map(sces, modelGeneVar)
combined.dec <- combineVar(decs)
hvg.var <- getTopHVGs(combined.dec, n=5000)

tes <- rowData(uncorrected0)[rowData(uncorrected0)$gene_biotype %in% "repeat_element",] |> rownames()
hvtes <- hvg.var[hvg.var %in% tes]
chosen <- hvg.var[!hvg.var %in% tes]
rowSubset(uncorrected0) <- chosen
metadata(uncorrected0)$highly.variable.tes <- hvtes
metadata(uncorrected0)$highly.variable.genes <- chosen

# ------------------------------------------------------------------------------
# UMAP + TSNE uncorrected - spoiler - def need batch correction
# ------------------------------------------------------------------------------

# PCA - without influence of mito genes or TEs
uncorrected <- runPCA(uncorrected0, subset_row=chosen, ncomponents=10)

g_pc_elbow_uncorrected <- tibble(percent.var = attr(reducedDim(uncorrected), "percentVar")) |>
  mutate(PC = row_number()) |>
  ggplot(aes(PC, percent.var)) +
  geom_point() +
  geom_vline(xintercept = ncol(reducedDim(uncorrected,"PCA")), color="red")

g_pca_uncorrected <- plotPCA(uncorrected, colour_by="batch")


set.seed(2)
uncorrected <- runUMAP(uncorrected, dimred="PCA",pca=5)
set.seed(2)
uncorrected <- runTSNE(uncorrected, dimred="PCA",pca=5)

g_umap_uncorrected <- plotReducedDim(uncorrected, dimred = "UMAP", colour_by = "batch")
g_tsne_uncorrected <- plotReducedDim(uncorrected, dimred = "TSNE", colour_by = "batch")

# ------------------------------------------------------------------------------
# integration
# https://bioconductor.org/books/3.17/OSCA.multisample/integrating-datasets.html#mnn-correction
# "...k... can be interpreted as the minimum anticipated frequency of any shared cell type or state in each batch."
# ------------------------------------------------------------------------------
set.seed(1000101001)
mnn.out <- fastMNN(uncorrected0, k=10,d = 50,
                   batch = uncorrected0$batch, 
                   subset.row=chosen, # no tes involved
                   correct.all = T, 
                   BSPARAM=BiocSingular::RandomParam(deferred = T))

corrected <- uncorrected0
reducedDims(corrected) <- reducedDims(mnn.out)
assay(corrected,"reconstructed") <- assay(mnn.out,"reconstructed")

# ------------------------------------------------------------------------------
# UMAP + TSNE corrected
# ------------------------------------------------------------------------------

g_pca_corrected <- plotReducedDim(corrected, dimred = "corrected",colour_by = "batch")

set.seed(2)
corrected <- runUMAP(corrected, dimred="corrected", pca=5)
set.seed(2)
corrected <- runTSNE(corrected, dimred="corrected", pca=5)

g_umap_corrected <- plotReducedDim(corrected, dimred = "UMAP", colour_by = "batch")
g_tsne_corrected <- plotReducedDim(corrected, dimred = "TSNE", colour_by = "batch")

# ------------------------------------------------------------------------------
# rerun clustering
# ------------------------------------------------------------------------------

colData(corrected) <- colData(corrected) |>
  as_tibble(rownames="ix") |>
  dplyr::rename_with(.fn=~{paste0("premerge_",.x)}, .cols = c("label","label2","celltype","celltype.garnett")) |>
  column_to_rownames("ix") |>
  DataFrame()

set.seed(20)
clust.sweep <- clusterSweep(reducedDim(corrected, "corrected"), 
                            NNGraphParam(), 
                            k=as.integer(c(35, 45, 55)),
                            cluster.fun=c("walktrap", "leiden"))

sweep.df <- as_tibble(clust.sweep$parameters,rownames="paramset")

sweep.df <- full_join(sweep.df, enframe(as.list(clust.sweep$clusters),name = "paramset", value = "labels"), by="paramset", relationship="one-to-one")

sweep.df <- sweep.df |>
  mutate(nclust = map_int(labels, ~length(unique(.x))))

sweep.df <- sweep.df |>
  mutate(mean.sil.width= map_dbl(labels,~mean(approxSilhouette(reducedDim(corrected,"corrected"),.x)$width)))

sweep.df <- sweep.df |>
  mutate(wcss= map_dbl(labels,~sum(clusterRMSD(reducedDim(corrected,"corrected"), .x, sum=TRUE), na.rm=TRUE)))

g_coarse_clustering_sweep <- sweep.df |>
  dplyr::select(-labels) |>
  pivot_longer(-c(paramset,k,cluster.fun),names_to = "metric", values_to = "value") |>
  ggplot(aes(k,value,color=cluster.fun)) +
  facet_wrap(~metric,scales = "free") +
  geom_line()


set.seed(20)
nn.clust <- clusterCells(corrected, use.dimred="corrected", full=TRUE,
                         BLUSPARAM=NNGraphParam(cluster.fun="leiden",k =45))

colLabels(corrected) <- nn.clust$clusters

plotReducedDim(corrected, dimred = "TSNE", colour_by = "label",text_by = "label") + paletteer::scale_color_paletteer_d("ggsci::default_igv")

sil.approx <- approxSilhouette(reducedDim(corrected, "corrected"), clusters=colLabels(corrected)) |>
  as_tibble(rownames = "Barcode")

sil.approx$closest <- factor(ifelse(sil.approx$width > 0, colLabels(corrected), sil.approx$other))
sil.approx$cluster <- colLabels(corrected)

g_silhouette <- ggplot(sil.approx, aes(x=cluster, y=width, colour=closest)) +
  geom_jitter()

# ------------------------------------------------------------------------------
# cluster naming 
# https://cole-trapnell-lab.github.io/garnett/docs_m3/#loading-your-data
# https://cole-trapnell-lab.github.io/garnett/docs_m3/#installing-garnett

mono_cds <- new_cell_data_set(logcounts(corrected),cell_metadata = colData(corrected),gene_metadata = rowData(corrected))

garnett_fl <- "data/green2018_garnett_adult_testis_classifier.txt"
garnett_fl <- snakemake@input$garnett_fl

marker_check <- check_markers(mono_cds, garnett_fl,
                              db=org.Mm.eg.db,
                              cds_gene_id_type = "ENSEMBL",
                              marker_file_gene_id_type = "SYMBOL")

g_plot_markers <- plot_markers(marker_check)

set.seed(20)
threads <- 8
threads <- snakemake@threads
garnett_classifier <- train_cell_classifier(cds = mono_cds,cores = threads,
                                            marker_file = garnett_fl,
                                            db=org.Mm.eg.db,
                                            cds_gene_id_type = "ENSEMBL",
                                            marker_file_gene_id_type = "SYMBOL")

mono_cds <- classify_cells(mono_cds, garnett_classifier,
                           db = org.Mm.eg.db,
                           cluster_extend = TRUE,
                           cds_gene_id_type = "ENSEMBL")


stopifnot(all(rownames(colData(mono_cds)) == rownames(colData(corrected))))

corrected$celltype <- colData(mono_cds)$cluster_ext_type

# ------------------------------------------------------------------------------
# cluster naming by vote - where some cells in a cluster have garnett labels that
# disagree with the predominant cluster type
colData(corrected) <- colData(corrected) |>
  as_tibble(rownames="cell") |>
  filter(celltype!="Unknown") |>
  dplyr::count(label, celltype) |>
  filter(!is.na(celltype)) |>
  group_by(label) |>
  slice_max(n,n=1,with_ties = F) |>
  ungroup() |> 
  dplyr::select(label, celltype) |>
  left_join(x=as_tibble(colData(corrected),rownames="cell"),y=_, by="label", suffix=c('.garnett',"")) |>
  column_to_rownames("cell") |>
  DataFrame()


corrected$label2 <- paste(corrected$label,corrected$celltype,sep="/")

# ------------------------------------------------------------------------------
# sanity checks
# ------------------------------------------------------------------------------

# plotReducedDim(corrected, dimred="UMAP",colour_by="celltype.garnett", text_by="label") +
#  facet_wrap(~colour_by)


# library(bluster)
# # better alignment (between within batch celltype assignment and post-merge) for germ cells, not
# # so good for somatic. likely due to far fewer high confidence somatic cells identified.
# cluster_nesting_9646 <- nestedClusters(ref=corrected[,corrected$batch=="batch9646"]$label2,
#                         alt=corrected[,corrected$batch=="batch9646"]$premerge_label2)
# 
# cluster_nesting_9646$alt.mapping |>
#   as_tibble(rownames="label2") |>
#   print(n=Inf)
# 
# 
# 
# lkup <- rowData(corrected) |> as_tibble(rownames("gene_id")) |> dplyr::select(gene_name, gene_id) |> deframe()
# plotReducedDim(corrected[,order(counts(corrected)[lkup["Tex14"],])],"UMAP", colour_by= lkup["Tex14"])
# plotReducedDim(corrected[,corrected$celltype %in% c("Spermatogonia","Elongating","RoundSpermatid","Spermatocyte")],"UMAP", colour_by= "celltype")
# plotReducedDim(corrected[,!corrected$celltype %in% c("Spermatogonia","Elongating","RoundSpermatid","Spermatocyte")],"UMAP", colour_by= "celltype")
# plotReducedDim(corrected,"UMAP", colour_by= "celltype")
# 
# osca suggests that this should be less than 10%
# metadata(mnn.out)$merge.info$lost.var
# 
# 
# # look at genes with high absolute batch effect
# common <- rownames(mnn.out)
# vars <- mnnDeltaVariance(sces$adult_9646_combined[common,], sces$adult_9682_combined[common,], 
#                          pairs=metadata(mnn.out)$merge.info$pairs)
# vars[order(vars$adjusted, decreasing=TRUE),]



write_rds(corrected,snakemake@output$rds)
write_rds(uncorrected, snakemake@output$uncorrected)
write_rds(garnett_classifier, snakemake@output$garnett_classifier)
write_rds(g_silhouette, snakemake@output$g_silhouette)
write_rds(g_coarse_clustering_sweep, snakemake@output$g_coarse_clustering_sweep)
write_rds(marker_check, snakemake@output$marker_check)
write_rds(mnn.out, snakemake@output$mnn_out)