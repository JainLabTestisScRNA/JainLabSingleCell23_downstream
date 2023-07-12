library(tidyverse)
library(scater)
library(scran)
library(scuttle)
library(ggpubr)
library(ggdensity)
library(patchwork)
library(paletteer)
library(tradeSeq)
library(pheatmap)
library(TSCAN)
library(slingshot)

sce_fl <- "results/single-cell-preproc/preprocessed/juvenile_13d_wt_null.cellranger.sce.germ_cell.rds"
sce_fl <- snakemake@input$sce
sce <- read_rds(sce_fl)

set.seed(2)
dec <- modelGeneVar(sce)
dec.tech <- modelGeneVarByPoisson(sce)
hvg.var <- getTopHVGs(dec, prop = 0.02)

tes <- rowData(sce)[rowData(sce)$gene_biotype %in% "repeat_element",] |> rownames()

hvtes <- hvg.var[hvg.var %in% tes]

#& 
chosen <- hvg.var[!hvg.var %in% tes & !hvg.var %in% subset(rowData(sce),grepl("^mt-",gene_name))$gene_name]

rowSubset(sce) <- chosen
metadata(sce)$highly.variable.tes <- hvtes
metadata(sce)$highly.variable.genes <- chosen

sce <- denoisePCA(sce, dec.tech, subset.row=chosen)

set.seed(2)
sce <- runUMAP(sce, n_neighbors=40, spread=3,  min_dist=2, n_trees=100, subset_row = chosen)

#plotReducedDim(sce, dimred = "UMAP", colour_by="sublabel2",text_by = "sublabel2")
#plotReducedDim(sce, dimred = "PCA", colour_by="sublabel2",text_by = "sublabel2")

# ------------------------------------------------------------------------------
# compute entropy
entropy <- perCellEntropy(sce,assay.type="logcounts")

sce$entropy <- entropy

g_entropy <- plotColData(sce, x = "sublabel2", y="entropy") +
  stat_summary() +
  theme(axis.text.x = element_text(angle=45, hjust=1))

#plotReducedDim(sce,dimred = "UMAP", colour_by = "entropy", text_by = "sublabel2")
# ------------------------------------------------------------------------------
# get pseudotime

# first step to ensure correct ordering
sce.sling0 <- slingshot(sce, reducedDim = "PCA", clusterLabels = factor(sce$sublabel2))
sce.sling0$slingshot <- NULL

reroot_df <- makePerCellDF(sce.sling0) |>
  as_tibble() |>
  group_by(label2) |>
  summarise(entropy=mean(entropy), ps = mean(slingPseudotime_1)) |>
  arrange(-ps)

root_label2 <- reroot_df[c(which.min(reroot_df$ps),which.max(reroot_df$ps)),] |>
  slice_max(entropy) |>
  pull(label2)

leaf_label2 <- reroot_df[c(which.max(reroot_df$ps),which.max(reroot_df$ps)),] |>
  slice_min(entropy) |>
  pull(label2)

sce.sling <- slingshot(sce, reducedDim = "PCA", clusterLabels = factor(sce$label2),start.clus = root_label2)

embedded <- embedCurves(sce.sling, "UMAP")
embedded <- slingCurves(embedded)[[1]] # only 1 path.
embedded <- data.frame(embedded$s[embedded$ord,])
embedded2 <- embedded[colnames(sce.sling),]

reducedDim(sce.sling, "sling") <- embedded2

saveRDS(sce.sling, snakemake@output$rds)

# colData(sce.sling)[,-which(colnames(colData(sce.sling)) == "slingshot")] |>
#   as_tibble() |>
# ggplot(aes(slingPseudotime_1, 1)) +
#   geom_jitter()
# 
# get_ps_umap <- \(x) plotUMAP(x, colour_by="slingPseudotime_1") +
#   geom_path(data=embedded[seq(nrow(embedded)-50),], aes(x=Dim.1, y=Dim.2), size=1.2, arrow = arrow(length=unit(0.30,"cm"), ends="last", type = "closed"))
# 
# get_ps_umap(sce.sling[,sce.sling$Sample==c("CMO306")]) +
# get_ps_umap(sce.sling[,sce.sling$Sample==c("CMO309")])

#ps_per_sample <- c("CMO306","CMO309")

#names(ps_per_sample) <- ps_per_sample

#map_df(ps_per_sample,~as_tibble(slingPseudotime(sce.sling[sce.sling$Sample==.x])),.id="Sample") |>
#  ggplot(aes(Lineage1,color=Sample)) +
#  geom_density()
#  facet_wrap(~Sample)

#ks.test(slingPseudotime(sce.sling)[colData(sce.sling)$Sample == "CMO306", 1],
#        slingPseudotime(sce.sling)[colData(sce.sling)$Sample == "CMO309", 1])

