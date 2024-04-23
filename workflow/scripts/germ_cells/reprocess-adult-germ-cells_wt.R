library(tidyverse)
library(DropletUtils)
library(scater)
library(scran)
library(scuttle)
library(bluster)
library(PCAtools)
library(batchelor)
library(TSCAN)

sce_fl <- "results/germ_cells/adult.sce.germ_cell.wt.Spermatogonia.rds"
sce_fl <-  snakemake@input$sce
sce <- read_rds(sce_fl)

k <- 2
k  <- snakemake@params$k


# ------------------------------------------------------------------------------
# variable feature selection
dec <- modelGeneVar(sce,block=sce$batch)
hvg.var <- getTopHVGs(dec, prop = 0.1)
  
tes <- rowData(sce)[rowData(sce)$gene_biotype %in% "repeat_element" ,] |> rownames()
  
hvtes <- hvg.var[hvg.var %in% tes]
  
chosen <- hvg.var[!hvg.var %in% tes & !hvg.var %in% subset(rowData(sce),grepl("^mt-",gene_name))$gene_name]
  
rowSubset(sce) <- chosen
metadata(sce)$highly.variable.tes <- hvtes
metadata(sce)$highly.variable.genes <- chosen
  
# ------------------------------------------------------------------------------
# PCA - without influence of TEs
sce <- fixedPCA(sce, subset.row=chosen,rank = 50)
  
percent.var <- attr(reducedDim(sce,"PCA"), "percentVar")

chosen.elbow <- findElbowPoint(percent.var)

sce <- fixedPCA(sce, subset.row=chosen,rank = chosen.elbow)
# ------------------------------------------------------------------------------
# integration
# https://bioconductor.org/books/3.17/OSCA.multisample/integrating-datasets.html#mnn-correction
# "...k... can be interpreted as the minimum anticipated frequency of any shared cell type or state in each batch."
# ------------------------------------------------------------------------------  
set.seed(2)
mnn.out <- fastMNN(sce, k=k,d = chosen.elbow,
                     batch = sce$batch, 
                     subset.row=chosen, # no tes involved
                     correct.all = T, 
                     BSPARAM=BiocSingular::RandomParam(deferred = T))
  
reducedDims(sce) <- reducedDims(mnn.out)
assay(sce,"reconstructed") <- assay(mnn.out,"reconstructed")
  
# ------------------------------------------------------------------------------
# recluster, shooting for k (smk param) clusters

set.seed(2)
nn.clust <- clusterCells(sce,use.dimred = "corrected", full=F,
                           BLUSPARAM=KmeansParam(k))
  
colLabels(sce) <- nn.clust |> as.character() |> paste(sce$celltype,sep="/")

plotReducedDim(sce,dimred = "corrected",ncomponents = c(1,2),colour_by = "label",swap_rownames = "gene_name")

write_rds(sce,snakemake@output$rds)



