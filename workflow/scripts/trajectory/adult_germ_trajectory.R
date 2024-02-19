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
library(batchelor)

sce_fl <- Sys.glob("results/single-cell-preproc/integrated/adult.cellranger.sce.integrated.rds")
sce_fl <- snakemake@input$sce

sce <- read_rds(sce_fl)
sce$Sample <- sce$genotype

sce <- sce[,sce$celltype%in% c("Elongating","RoundSpermatid","Spermatocyte","Spermatogonia")]

# ------------------------------------------------------------------------------
# lightweight integration for germ cells only
# slingshot consistently errors w/ the fastMNN 'reconstructed' passed as the dimred
# so do multibatchPCA instead
# ------------------------------------------------------------------------------
assay(sce,"reconstructed") <- NULL
reducedDims(sce) <- NULL

# get hvgs - recompute for germ only for higher granularity
dec <- modelGeneVar(sce,block=sce$batch)
hvg.var <- getTopHVGs(dec, n=5000)

tes <- rowData(sce)[rowData(sce)$gene_biotype %in% "repeat_element",] |> rownames()
hvtes <- hvg.var[hvg.var %in% tes]
chosen <- hvg.var[!hvg.var %in% tes]
rowSubset(sce) <- chosen
metadata(sce)$highly.variable.tes <- hvtes
metadata(sce)$highly.variable.genes <- chosen

stopifnot(all(rownames(sce)[rowSubset(sce)] %in% chosen))

pcx <- multiBatchPCA(sce,batch=sce$batch,subset.row = rowSubset(sce))

reducedDim(sce,"PCA") <- rbind(pcx$batch9682,pcx$batch9646)[colnames(sce),]

set.seed(2)
sce <- runUMAP(sce, n_neighbors=25, spread=2,  min_dist=1, n_trees=75, n_dimred=5,dimred="PCA")

plotReducedDim(sce, dimred = "UMAP", colour_by="batch",text_by = "label2",swap_rownames="gene_name")

# ------------------------------------------------------------------------------
# compute entropy
entropy <- perCellEntropy(sce,assay.type="logcounts",)

sce$entropy <- entropy

g_entropy <- plotColData(sce, x = "label2", y="entropy",other_fields = c("batch","genotype")) +
  stat_summary(position = position_dodge(width=0.5)) + 
  aes(color=genotype) +
  theme(axis.text.x = element_text(angle=45, hjust=1))


# ------------------------------------------------------------------------------
# get pseudotime

# first step to ensure correct ordering
#plotReducedDim(sce,dimred = "UMAP", colour_by = "entropy", text_by = "label2")

# based on plot above, specify first/last cluster
#https://github.com/kstreet13/slingshot/issues/37
sce.sling <- slingshot(sce, reducedDim = "PCA", clusterLabels = sce$celltype, start.clus="Spermatogonia",end.clus="Elongating")

embedded <- embedCurves(sce.sling, "UMAP")
embedded <- slingCurves(embedded)[[1]] # only 1 path.
embedded <- data.frame(embedded$s[embedded$ord,])
embedded2 <- embedded[colnames(sce.sling),]

reducedDim(sce.sling, "sling") <- embedded2

saveRDS(sce.sling, snakemake@output$rds)

#colData(sce.sling)[,-which(colnames(colData(sce.sling)) == "slingshot")] |>
#  as_tibble() |>
#  ggplot(aes(slingPseudotime_1, 1)) +
#  geom_jitter() +
#  facet_wrap(~genotype)

#get_ps_umap <- \(x) plotUMAP(x, colour_by="slingPseudotime_1") +
#  geom_path(data=embedded[seq(nrow(embedded)-50),], aes(x=Dim.1, y=Dim.2), size=1.2, arrow = arrow(length=unit(0.30,"cm"), ends="last", type = "closed")) +
#  labs(title=unique(x$genotype))

#get_ps_umap(sce.sling[,sce.sling$genotype==c("WT")]) +
#get_ps_umap(sce.sling[,sce.sling$genotype==c("MUT")])

#ps_per_sample <- c("MUT","WT")

#names(ps_per_sample) <- ps_per_sample

#map_df(ps_per_sample,~as_tibble(slingPseudotime(sce.sling[,sce.sling$genotype==.x])),.id="genotype") |>
#  ggplot(aes(Lineage1,fill=genotype)) +
#  geom_histogram() +
#  facet_wrap(~genotype)
