library(tidyverse)
library(scater)
library(scran)
library(scuttle)
library(bluster)

sce <- read_rds("results/single-cell-preproc/integrated/adult.cellranger.sce.passing_cells.rds")

lkup <- rowData(sce) |> as_tibble(rownames("gene_id")) |> dplyr::select(gene_name, gene_id) |> deframe()
goi <- lkup["Ddx4"]

plotReducedDim(sce[,order(counts(sce)[lkup["Ddx4"],])],"UMAP", colour_by= lkup["Ddx4"])
plotReducedDim(sce[,order(counts(sce)[lkup["Tex14"],])],"UMAP", colour_by= lkup["Tex14"])
plotReducedDim(sce[,order(counts(sce)[lkup["Sycp3"],])],"UMAP", colour_by= lkup["Sycp3"])
plotReducedDim(sce[,order(counts(sce)[lkup["Sycp3"],])],"UMAP", colour_by= lkup["Dnmt3L"])
plotReducedDim(sce[,order(counts(sce)[lkup["Tex101"],])],"UMAP", colour_by= lkup["Tex101"])


plotExpression(sce, features=gois, x=label)


plotReducedDim(sce[,order(counts(sce)[lkup["Tex14"],])],"UMAP", colour_by= lkup["Tex14"],other_fields = "batch") + facet_wrap(~batch)

plotReducedDim(sce[,order(counts(sce)[lkup["Sycp3"],])],"UMAP", colour_by= lkup["Sycp3"],other_fields = "batch") + facet_wrap(~batch)

plotReducedDim(sce[,order(counts(sce)[lkup["Ddx4"],])],"UMAP", colour_by= lkup["Ddx4"],other_fields = "batch") + facet_wrap(~batch)

plotReducedDim(sce[,order(sce$sum)],"UMAP", colour_by= "sum",other_fields = "batch") + facet_wrap(~batch)

plotReducedDim(sce[,sce$sum > 2000],"UMAP", colour_by= "sum",other_fields = "batch") + facet_wrap(~batch)

# ------------------------------------------------------------------------------
# coarse clustering
# set.seed(20)
# clust.sweep <- clusterSweep(reducedDim(sce, "corrected"), 
#                             NNGraphParam(), 
#                             k=as.integer(c(15, 25, 35, 45, 55, 65, 75, 85)),
#                             cluster.fun=c("louvain", "walktrap", "infomap","leiden"),
#                             BPPARAM=BiocParallel::MulticoreParam(8))
# 
# sweep.df <- as_tibble(clust.sweep$parameters,rownames="paramset")
# 
# sweep.df <- full_join(sweep.df, enframe(as.list(clust.sweep$clusters),name = "paramset", value = "labels"), by="paramset", relationship="one-to-one")
# 
# sweep.df <- sweep.df |>
#   mutate(nclust = map_int(labels, ~length(unique(.x))))
# 
# sweep.df <- sweep.df |>
#   mutate(mean.sil.width= map_dbl(labels,~mean(approxSilhouette(reducedDim(sce,"PCA"),.x)$width)))
# 
# sweep.df <- sweep.df |>
#   mutate(wcss= map_dbl(labels,~sum(clusterRMSD(reducedDim(sce), .x, sum=TRUE), na.rm=TRUE)))
# 
# g_coarse_clustering_sweep <- sweep.df |>
#   dplyr::select(-labels) |>
#   pivot_longer(-c(paramset,k,cluster.fun),names_to = "metric", values_to = "value") |>
#   ggplot(aes(k,value,color=cluster.fun)) +
#   facet_wrap(~metric,scales = "free") +
#   geom_line()


set.seed(20)
nn.clust <- clusterCells(sce, use.dimred="corrected", full=TRUE, assay.type=NULL,
                         BLUSPARAM=NNGraphParam(cluster.fun="leiden",k =125))

colLabels(sce) <- nn.clust$clusters

plotReducedDim(sce, dimred = "UMAP", colour_by = "label",text_by = "label") + paletteer::scale_color_paletteer_d("ggsci::default_igv")
#plotReducedDim(sce, dimred = "UMAP", colour_by = "x",text_by = "x") + paletteer::scale_color_paletteer_d("ggsci::default_igv")


goi <- c("Tex14", "Sycp3","Ddx4") #|> map_chr(~{lkup[.x]})
plotExpression(sce,x = "label",features = goi)


sil.approx <- approxSilhouette(reducedDim(sce, "PCA"), clusters=colLabels(sce)) |>
  as_tibble(rownames = "Barcode")

sil.approx$closest <- factor(ifelse(sil.approx$width > 0, colLabels(sce), sil.approx$other))
sil.approx$cluster <- colLabels(sce)

g_silhouette <- ggplot(sil.approx, aes(x=cluster, y=width, colour=closest)) +
  geom_jitter()