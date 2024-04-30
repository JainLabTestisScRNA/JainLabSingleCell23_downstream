Sys.setenv(R_PROFILE=".Rprofile")
source(Sys.getenv("R_PROFILE"))

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

fl <- "results/integration/adult.sce.integrated.rds"
fl <- snakemake@input$sce

corrected <- read_rds(fl)

# ------------------------------------------------------------------------------
# rerun clustering
# ------------------------------------------------------------------------------
# set.seed(20)
# clust.sweep <- clusterSweep(reducedDim(corrected, "corrected"), 
#                             SNNGraphParam(), 
#                             k=as.integer(c(25,35, 45, 55,65)),
#                             cluster.fun=c("louvain"),type="jaccard")
# 
# sweep.df <- as_tibble(clust.sweep$parameters,rownames="paramset")
# 
# sweep.df <- full_join(sweep.df, enframe(as.list(clust.sweep$clusters),name = "paramset", value = "labels"), by="paramset", relationship="one-to-one")
# 
# sweep.df <- sweep.df |>
#   mutate(nclust = map_int(labels, ~length(unique(.x))))
# 
# sweep.df <- sweep.df |>
#   mutate(mean.sil.width= map_dbl(labels,~mean(approxSilhouette(reducedDim(corrected,"corrected"),.x)$width)))
# 
# sweep.df <- sweep.df |>
#   mutate(wcss= map_dbl(labels,~sum(clusterRMSD(reducedDim(corrected,"corrected"), .x, sum=TRUE), na.rm=TRUE)))
# 
# g_coarse_clustering_sweep <- sweep.df |>
#   dplyr::select(-labels) |>
#   pivot_longer(-c(paramset,k,cluster.fun),names_to = "metric", values_to = "value") |>
#   ggplot(aes(k,value,color=cluster.fun)) +
#   facet_wrap(~metric,scales = "free") +
#   geom_line()


set.seed(20)
nn.clust <- clusterCells(corrected, use.dimred="corrected", full=TRUE,
                         BLUSPARAM=NNGraphParam(cluster.fun="leiden",k =55))

colLabels(corrected) <- nn.clust$clusters

plotReducedDim(corrected, dimred = "UMAP", colour_by = "label",text_by = "label") + paletteer::scale_color_paletteer_d("ggsci::default_igv")

sil.approx <- approxSilhouette(reducedDim(corrected, "corrected"), clusters=colLabels(corrected)) |>
  as_tibble(rownames = "Barcode")

sil.approx$closest <- factor(ifelse(sil.approx$width > 0, colLabels(corrected), sil.approx$other))
sil.approx$cluster <- colLabels(corrected)

g_silhouette <- ggplot(sil.approx, aes(x=cluster, y=width, colour=closest)) +
  geom_jitter()


write_rds(corrected, snakemake@output$rds)
write_rds(g_silhouette, snakemake@output$g_silhouette)
#write_rds(g_coarse_clustering_sweep, snakemake@output$g_coarse_clustering_sweep)