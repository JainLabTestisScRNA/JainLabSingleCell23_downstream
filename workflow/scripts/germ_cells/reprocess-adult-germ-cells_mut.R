library(tidyverse)
library(DropletUtils)
library(scater)
library(scran)
library(scuttle)
library(bluster)
library(PCAtools)

sce_fl <- "results/germ_cells/adult.sce.integrated.clustered.celltypes.germ_cell.mut.rds"
sce_fl <-  snakemake@input$sce
sce <- read_rds(sce_fl)
assay(sce,"reconstructed") <- NULL
reducedDims(sce) <- NULL
sce$label2 <- NULL
sce$celltype.garnett <- NULL

# see methods of green et al 'Focused analysis for germ cells'
sce <- sce[,sce$nexprs >= 1000]

# ------------------------------------------------------------------------------
# variable feature selection
dec <- modelGeneVar(sce,block=sce$batch)
hvg.var <- getTopHVGs(dec, prop=0.025)

tes <- rowData(sce)[rowData(sce)$gene_biotype %in% "repeat_element" ,] |> rownames()

hvtes <- hvg.var[hvg.var %in% tes]

chosen <- hvg.var[!hvg.var %in% tes & !hvg.var %in% subset(rowData(sce),grepl("^mt-",gene_name))$gene_name]

rowSubset(sce) <- chosen
metadata(sce)$highly.variable.tes <- hvtes
metadata(sce)$highly.variable.genes <- chosen

# ------------------------------------------------------------------------------
# PCA - without influence of TEs
sce <- fixedPCA(sce, subset.row=chosen,rank = 50)

# ------------------------------------------------------------------------------
# integration
# https://bioconductor.org/books/3.17/OSCA.multisample/integrating-datasets.html#mnn-correction
# "...k... can be interpreted as the minimum anticipated frequency of any shared cell type or state in each batch."
# ------------------------------------------------------------------------------
library(batchelor)

set.seed(2)
mnn.out <- fastMNN(sce, k=12,d = 50,
                   batch = sce$batch, 
                   subset.row=chosen, # no tes involved
                   correct.all = T, 
                   BSPARAM=BiocSingular::RandomParam(deferred = T))

reducedDims(sce) <- reducedDims(mnn.out)
assay(sce,"reconstructed") <- assay(mnn.out,"reconstructed")

# ------------------------------------------------------------------------------
# recluster, shooting for ~12 clusts
# want to be able to easily sanity check against Green et al. 2018

set.seed(2)
nn.clust <- clusterCells(sce,use.dimred = "corrected", full=F,
                         BLUSPARAM=SNNGraphParam(k = 48,type = "jaccard",cluster.fun = "louvain"))

colLabels(sce) <- nn.clust |> as.character()



plotReducedDim(sce[,order(sce$label)], dimred = "corrected", ncomponents = c(1,2),
               swap_rownames = "gene_name",other_fields = "genotype",
               colour_by = "label",
               text_by = "label") #+ facet_wrap(~genotype)

plotReducedDim(sce[,order(sce$label)], dimred = "corrected", ncomponents = c(1,2),
               swap_rownames = "gene_name",other_fields = "genotype",
               colour_by = "celltype",
               text_by = "label") #+ facet_wrap(~genotype)

# ------------------------------------------------------------------------------
# call cell types again, but based on Green 2018's 12 GC clusts
# https://cole-trapnell-lab.github.io/garnett/docs_m3/#loading-your-data
# https://cole-trapnell-lab.github.io/garnett/docs_m3/#installing-garnett

#library(garnett)
#library(org.Mm.eg.db)

#garnett_fl <- "data/green2018_garnett_adult_germcell_classifier.txt"
#garnett_fl <- snakemake@input$garnett_fl

# see docs for garnett's classify_cells.R
# garnett cluster column used for assignment by vote that I had previously done
#colData(sce)$garnett_cluster <- sce$label |> as.character()
#mono_cds <- new_cell_data_set(logcounts(sce),cell_metadata = colData(sce),gene_metadata = rowData(sce))

#marker_check <- check_markers(mono_cds, garnett_fl,
#                              db=org.Mm.eg.db,
#                              cds_gene_id_type = "ENSEMBL",
#                              marker_file_gene_id_type = "SYMBOL")

#g_plot_markers <- plot_markers(marker_check)

#set.seed(20)
#threads <- 4
#threads <- snakemake@threads
#garnett_classifier <- train_cell_classifier(cds = mono_cds,cores = threads,
#                                            marker_file = garnett_fl,
#                                           db=org.Mm.eg.db,
#                                           cds_gene_id_type = "ENSEMBL",
#                                             marker_file_gene_id_type = "SYMBOL")

# my instinct is that this is attempting to work on some very fine-grained clusters
# with substantial similarity amongst them. I tuned the params
# to tolerate much more uncertainty within clusters for cluster extension:
# because we get many unknowns and many clusters with no majority assignment
#mono_cds <- classify_cells(mono_cds, garnett_classifier,
#                           cluster_extend_max_frac_unknown = 0.97,
#                           cluster_extend_max_frac_incorrect=0.8,
#                           rank_prob_ratio=1.5,
#                           db = org.Mm.eg.db,
#                           cluster_extend = TRUE,
#                          cds_gene_id_type = "ENSEMBL")


#stopifnot(all(rownames(colData(mono_cds)) == rownames(colData(sce))))

#sce$celltype <- colData(mono_cds)$cluster_ext_type


order_df <- tibble(label=(as.character(c(3,2,1,4,6,9,8,5,7)))) |>
  mutate(ordered_label = paste0("cl",row_number()))

colData(sce) <- colData(sce) |> as_tibble(rownames="cell") |>
  left_join(order_df,by=c("label")) |>
  column_to_rownames("cell") |>
  DataFrame()

plotReducedDim(sce,"corrected",ncomponents = c(1,2),colour_by = "ordered_label",text_by = "ordered_label")
sce$clustering_label = sce$label
sce$label <- sce$ordered_label


write_rds(sce,snakemake@output$rds)



