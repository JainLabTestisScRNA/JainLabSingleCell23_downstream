library(tidyverse)
library(scater)
library(scran)
library(scuttle)
library(ggpubr)
library(ggdensity)
library(patchwork)
library(paletteer)
library(dtwclust)
library(tradeSeq)
library(pheatmap)
library(TSCAN)
library(slingshot)

library(clusterExperiment)

sce_fl <- "results/tradeseq/juvenile_13d_wt_null.cellranger.sce.germ_cell.tradeseq.rds"
res_fl <- "results/tradeseq/juvenile_13d_wt_null.cellranger.germ_cell.deg.rds"

sce_fl <- snakemake@input$sce
res_fl <- snakemake@input$deg

sce <- read_rds(sce_fl)
conditionRes <- read_rds(res_fl)

conditionGenes <- conditionRes$feature[conditionRes$padj <= 0.1]
conditionGenes <- conditionGenes[!is.na(conditionGenes)]


id_lkup <- rowData(sce)[,c("gene_id", "gene_name")] |>
  as_tibble() |>
  dplyr::rename(feature="gene_id")

# get predicted (non transformed) values and scale them
yhatSmooth <- predictSmooth(sce, gene = conditionGenes, nPoints = 100, tidy = T) |> 
  as_tibble() |> 
  group_by(lineage,gene, condition) |>
  mutate(yhatScaled = scale(yhat)[,1])

yhatSmoothScaledMat<- predictSmooth(sce, gene = conditionGenes, nPoints = 100, tidy = F) |> t() |> scale() |> t()

res <- yhatSmooth |>
  arrange(time) |>
  group_by(gene,lineage,condition) |>
  mutate(pseudotime.bin = row_number()) |>
  ungroup() |>
  relocate(pseudotime.bin, .after="time") |>
  dplyr::rename(feature = "gene")

res <- left_join(res, id_lkup, by="feature") |> relocate(gene_name, .after = "feature")

res <- res |>
  dplyr::select(lineage, feature,yhat, condition, pseudotime.bin) |>
  pivot_wider(names_from = "condition", values_from = "yhat") |>
  mutate(lfcMutant.WT = log2((MUT+1)/(WT+1))) |>
  dplyr::select(lineage, feature, pseudotime.bin, lfcMutant.WT) |>
  left_join(res,y=_, by=c("lineage","feature","pseudotime.bin"))


# ------------------------------------------------------------------------------
# code interlude - choose k and do kmeans

get_silh <- \(x, y) bluster::approxSilhouette(x, y$cluster[rownames(x)]) |> as_tibble(rownames="feature")
get_rmsd <- \(x, y) sum(bluster::clusterRMSD(x, y$cluster[rownames(x)], sum=T), na.rm = T)


estimate_k <- function(mat,ks = 3:20) {
  set_names(ks) |>
    map_df(~{
      km <- kmeans(mat, .x)
      set.seed(2)
      sil <- get_silh(mat, km)
      mutate(sil, rmsd = get_rmsd(mat, km))
    }, .id="k")
}

estimated_k_tbl <- estimate_k(yhatSmoothScaledMat) |>
  group_by(k) |>
  summarise(width = mean(width),rmsd = unique(rmsd)) |>
  mutate(k=fct_reorder(factor(k),as.integer(k)))

k_choice <- estimated_k_tbl |>
  mutate(width_rank = dense_rank(-width), rmsd_rank = dense_rank(rmsd)) |>
  group_by(k) |>
  summarise(score = width_rank + rmsd_rank) |>
  slice_min(score) |>
  slice_min(k) |>
  pull(k) |>
  as.character() |>
  as.integer()

g_k_choice <- estimated_k_tbl |>
  pivot_longer(-k) |>
  ggplot(aes(k, value)) +
  geom_point() +
  facet_wrap(~name,scales = "free")

simple_clust <- yhatSmoothScaledMat |>
  kmeans(k_choice)

# ------------------------------------------------------------------------------
# back to analysis


res <- enframe(simple_clust$cluster, name = "feature", value = "cluster") |>
  mutate(cluster=as.character(cluster)) |>
  left_join(res,y=_, by="feature")

res <-  res |> 
  mutate(gene_name = fct_reorder(gene_name, lfcMutant.WT))

clusters_df <- res |>
  dplyr::select(feature, gene_name, cluster) |>
  distinct()

saveRDS(res, snakemake@output$smoothed_cluster_timecourse_rds)
saveRDS(g_k_choice, snakemake@output$g_k_choice)
write_tsv(clusters_df, snakemake@output$gene_clusters_tsv)




