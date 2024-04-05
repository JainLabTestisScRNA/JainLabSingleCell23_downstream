library(tidyverse)
library(DropletUtils)
library(scater)
library(scran)
library(scuttle)
library(bluster)
library(PCAtools)
library(batchelor)

sce_fl <- "results/germ_cells/adult.sce.integrated.clustered.celltypes.germ_cell.rds"
sce_fl <-  snakemake@input$sce
sce <- read_rds(sce_fl)
assay(sce,"reconstructed") <- NULL
reducedDims(sce) <- NULL
sce$label2 <- NULL
sce$celltype.garnett <- NULL

# see methods of green et al 'Focused analysis for germ cells'
sce <- sce[,sce$nexprs >= 1000]
rowSubset(sce) <- NULL

sce.spg <- sce[,sce$celltype == "Spermatogonia"]
sce.spc <- sce[,sce$celltype == "Spermatocyte"]
sce.rspt <- sce[,sce$celltype == "RoundSpermatid"]
sce.espt <- sce[,sce$celltype == "Elongating"]



my_subcluster <- function(sce, k,ct) {
  # ------------------------------------------------------------------------------
  # variable feature selection
  dec <- modelGeneVar(sce,block=sce$batch)
  hvg.var <- getTopHVGs(dec, prop=0.5)
  
  tes <- rowData(sce)[rowData(sce)$gene_biotype %in% "repeat_element" ,] |> rownames()
  
  hvtes <- hvg.var[hvg.var %in% tes]
  
  chosen <- hvg.var[!hvg.var %in% tes & !hvg.var %in% subset(rowData(sce),grepl("^mt-",gene_name))$gene_name]
  
  rowSubset(sce) <- chosen
  metadata(sce)$highly.variable.tes <- hvtes
  metadata(sce)$highly.variable.genes <- chosen
  
  # ------------------------------------------------------------------------------
  # PCA - without influence of TEs
  sce <- fixedPCA(sce, subset.row=chosen,rank = 24)

  
  # ------------------------------------------------------------------------------
  # integration
  # https://bioconductor.org/books/3.17/OSCA.multisample/integrating-datasets.html#mnn-correction
  # "...k... can be interpreted as the minimum anticipated frequency of any shared cell type or state in each batch."
  # ------------------------------------------------------------------------------  
  set.seed(2)
  mnn.out <- fastMNN(sce, k=12,d = 24,
                     batch = sce$batch, 
                     subset.row=chosen, # no tes involved
                     correct.all = T, 
                     BSPARAM=BiocSingular::RandomParam(deferred = T))
  
  reducedDims(sce) <- reducedDims(mnn.out)
  assay(sce,"reconstructed") <- assay(mnn.out,"reconstructed")
  
  # ------------------------------------------------------------------------------
  # recluster, shooting for ~12 clusts
  # want to be able to easily sanity check against Green et al. 2018
  
  colLabels(sce) <- sce$celltype
  
  set.seed(2)
  nn.clust <- clusterCells(sce,use.dimred = "corrected", full=F,
                           BLUSPARAM=SNNGraphParam(type="jaccard",cluster.fun="louvain",k=k))
  
  colLabels(sce) <- nn.clust |> as.character() |> paste(ct,sep="/")
  
  rowSubset(sce) <- NULL
  assay(sce,"reconstructed") <- NULL
  reducedDims(sce) <- NULL
  
  sce
}

sce.spg$label <- "1/Spermatogonia"

# special subclust of spermatocytes
sce.spc0 <- sce.spc |> my_subcluster(12*100,"Spermatocyte")
sce.spc1 <- sce.spc0[,sce.spc0$label == "1/Spermatocyte"] |> my_subcluster(144,"SpcA")
sce.spc2 <- sce.spc0[,sce.spc0$label == "2/Spermatocyte"] |> my_subcluster(96,"SpcB")

sce.rspt <- sce.rspt |> my_subcluster(144*4,"RoundSpermatid")
sce.espt$label <- "1/Elongating"

rns <- list(sce.spg,sce.spc1,sce.spc2,sce.rspt,sce.espt) |>
  map(rownames)

# sanity check on rows to ensure suitable for cbind -  this is pure paranoia
stopifnot(all(all(rns[[1]] == rns[[2]]) &
all(rns[[1]] == rns[[3]]) &
all(rns[[1]] == rns[[4]]) &
all(rns[[1]] == rns[[5]])))

sce <- cbind(sce.spg,sce.spc1,sce.spc2,sce.rspt,sce.espt)

# variable feature selection
dec <- modelGeneVar(sce,block=sce$batch)
hvg.var <- getTopHVGs(dec, prop=0.1)

tes <- rowData(sce)[rowData(sce)$gene_biotype %in% "repeat_element" ,] |> rownames()

hvtes <- hvg.var[hvg.var %in% tes]

chosen <- hvg.var[!hvg.var %in% tes & !hvg.var %in% subset(rowData(sce),grepl("^mt-",gene_name))$gene_name]

rowSubset(sce) <- chosen
metadata(sce)$highly.variable.tes <- hvtes
metadata(sce)$highly.variable.genes <- chosen

# ------------------------------------------------------------------------------
# PCA - without influence of TEs
sce <- fixedPCA(sce, subset.row=chosen,rank = 24)


# ------------------------------------------------------------------------------
# integration
# https://bioconductor.org/books/3.17/OSCA.multisample/integrating-datasets.html#mnn-correction
# "...k... can be interpreted as the minimum anticipated frequency of any shared cell type or state in each batch."
# ------------------------------------------------------------------------------  
set.seed(2)
mnn.out <- fastMNN(sce, k=12,d = 50,
                   batch = sce$batch, 
                   subset.row=chosen, # no tes involved
                   correct.all = T, 
                   BSPARAM=BiocSingular::RandomParam(deferred = T))

reducedDims(sce) <- reducedDims(mnn.out)
assay(sce,"reconstructed") <- assay(mnn.out,"reconstructed")


plotReducedDim(sce[,order(sce$label)], dimred = "corrected", ncomponents = c(1,3),
              swap_rownames = "gene_name",other_fields = "genotype",
               colour_by = "label",
               text_by = "label") #+ facet_wrap(~genotype)


order_df <- tibble(label=c(
  "1/Spermatogonia",
  "4/SpcA",
  "1/SpcA",
  "3/SpcA",
  "4/SpcB",
  "3/SpcB",
  "1/SpcB",
  "2/SpcB",
  "2/SpcA",
  "3/RoundSpermatid",
  "2/RoundSpermatid",
  "1/RoundSpermatid",
  "1/Elongating"
)) |>
  mutate(ordered_label = paste0("cl",row_number()))

colData(sce) <- colData(sce) |> as_tibble(rownames="cell") |>
  left_join(order_df,by=c("label")) |>
  column_to_rownames("cell") |>
  DataFrame()

plotReducedDim(sce,"corrected",ncomponents = c(1,2),colour_by = "ordered_label",text_by = "ordered_label")
sce$clustering_label = sce$label
sce$label <- sce$ordered_label


write_rds(sce,snakemake@output$rds)



