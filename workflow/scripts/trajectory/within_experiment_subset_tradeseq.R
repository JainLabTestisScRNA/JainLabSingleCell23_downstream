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

threads <- ifelse(exists("snakemake"), snakemake@threads, 8)
knots <- ifelse(exists("snakemake"), snakemake@params$knots, 5)

sce_fl <- "results/trajectory/juvenile_13d_wt_null.cellranger.sce.germ_cell.trajectory.rds"
sce_fl <- snakemake@input$sce

sce <- readRDS(sce_fl)

stopifnot(all(metadata(sce)$highly.variable.genes %in% rownames(sce)[rowSubset(sce)]))

# running this on everything takes a very long time, so we'll
# make a very light cut.
dec <- modelGeneVar(sce)
hvg.var <- getTopHVGs(dec, prop = 0.5)

mito <- rownames(sce)[grep("chrM",as.character(seqnames(sce)))]

hvf <- hvg.var[!hvg.var %in% mito]

hvf_ix <- which(rownames(sce) %in% hvf)

# the vignette (https://kstreet13.github.io/bioc2020trajectories/articles/workshopTrajectories.html)
#set.seed(3)
#icMat <- evaluateK(counts = as.matrix(assays(sce)$counts),
#                   pseudotime = sce$slingPseudotime_1,
#                   cellWeights = pathStat(sce$slingshot,i="weights"),
#                  conditions = factor(sce$Sample),
#                  nGenes = 300,
#                   k = 3:8)

sce <- fitGAM(counts = sce,
              #pseudotime = sce$slingPseudotime_1,
              #cellWeights = pathStat(sce$slingshot,i="weights"), 
              conditions = factor(sce$Sample), 
              nknots = knots, 
              BPPARAM = BiocParallel::MulticoreParam(threads),
              genes = hvf_ix)

write_rds(sce, snakemake@output$rds)

