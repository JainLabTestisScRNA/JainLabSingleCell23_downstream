library(tidyverse)
library(DropletUtils)
library(scater)
library(scran)
library(scuttle)
library(ggdensity)
library(scDblFinder)
library(flexmix)
library(celda)
library(bluster)

# ------------------------------------------------------------------------------
# files

sce_fl <- "results/single-cell-preproc/import/juvenile_13d_wt_null.cellranger.sce.rds"
sce_fl <- snakemake@input$sce

raw_fl <- "results/single-cell-preproc/raw/juvenile_13d_wt_null.cellranger.sce.rds"
raw_fl <- snakemake@input$raw

celltype_fl <- "data/celltype_classification_per_barcode.csv"
celltype_fl <- snakemake@input$celltype_calls

cmo_fl <- "data/cellranger/custom/assignment_confidence_table.csv"
cmo_fl <- snakemake@input$sample_assignments
# ------------------------------------------------------------------------------
# import sce objs + add known cell types

sce <- read_rds(sce_fl)
raw <- read_rds(raw_fl) # for ambient rna removal step

celltype <- read_csv(celltype_fl)

colData(sce) <- as_tibble(colData(sce),rownames="cell") |>
  left_join(celltype, by=c("Barcode"="...1"), suffix=c("",".umich")) |>
  column_to_rownames("cell") |>
  DataFrame()

stopifnot(length(intersect(sce[,sce$Assignment=="CMO306"]$Barcode,
                           sce[,sce$Assignment=="CMO309"]$Barcode)) == 0)

#-------------------------------------------------------------------------------
# ambient rna removal

# use this to identify blank cells
cmo <- read_csv(cmo_fl)
non_empty_cmo <- cmo |> filter(!Assignment %in% c("Blank"))

set.seed(2)
dcx.sce <- decontX(sce, background=raw[,!raw$Barcode %in% non_empty_cmo$Barcode], seed=12345)

counts(sce) <- round(decontXcounts(dcx.sce))

g_dcx <- celda::plotDecontXContamination(dcx.sce)

# ------------------------------------------------------------------------------
# refilter empty drops after ambient rna removal

#sum(colSums(counts(sce)) > attr(discard.sum.madlower,"thresholds")[["lower"]])
#hist(colSums(counts(sce)))
discard.sum.madlower <- isOutlier(colSums(counts(sce)), type="lower", nmads=0.5, log = F)


set.seed(100)
e.out <- emptyDrops(counts(sce),lower = attr(discard.sum.madlower,"thresholds")[["lower"]])

# <0.0001 == a real cell
summary(e.out$FDR <= 0.001)

# <0.0001 == a real cell
metadata(sce)$empty_after_decontx <- e.out |>
  as_tibble() |>
  dplyr::count(realcell = FDR <= 0.001) |>
  filter(!realcell | is.na(realcell)) |>
  pull(n) |>
  sum()

sce <- sce[,which(e.out$FDR <= 0.001)]

# ------------------------------------------------------------------------------
# doublet removal

# get logcounts
sce <- logNormCounts(sce,assay.type="counts")

# find doublets
set.seed(2)
scdbf <- scDblFinder(sce)

# get doublet finder data as tbl
dat_doublet <- colData(scdbf) %>% as_tibble(rownames = "cell")

metadata(sce)$doublets_removed <- dat_doublet |>
  filter(scDblFinder.class == "doublet") |>
  nrow()

# perform filt
sce <- sce[,scdbf$scDblFinder.class == "singlet"]

# ------------------------------------------------------------------------------
# basic dead/dying qc
is.mito <- grep("chrM",as.character(seqnames(sce)))

is.te <- rowData(sce)$gene_biotype %in% c("repeat_element") & !rowData(sce)$gene_biotype %in% c("7SLRNA")

sce <- addPerCellQCMetrics(sce,subsets=(list(MT=is.mito,
                                             TE=is.te,
                                             `7SLRNA`=names(sce)[str_detect(names(sce),"7SL")])))

discard.mito.madTop <- isOutlier(sce$subsets_MT_percent, type="higher", nmads=2)
discard.mito.mad1 <- isOutlier(sce$subsets_MT_percent, type="higher", nmads=1)

g_mito_cuts <- tibble(total.decont.umis = sce$sum,
       mito.percent = sce$subsets_MT_percent) |>
  ggplot(aes(total.decont.umis,mito.percent)) +
  geom_hline(yintercept = attr(discard.mito.madTop,"thresholds")["higher"], color="red") +
  geom_hline(yintercept = attr(discard.mito.mad1,"thresholds")["higher"], color="blue") +
  geom_point() +
  scale_x_log10()

sce$mito_warning <- as.vector(discard.mito.mad1)
sce <- sce[,!discard.mito.madTop]

metadata(sce)$n_mito_discarded <- sum(discard.mito.madTop)
metadata(sce)$n_mito_warning <- sum(sce$mito_warning)

sce_mito_warning <- sce[,sce$mito_warning]
sce <- sce[,!sce$mito_warning]

# ------------------------------------------------------------------------------
# remove low dimension cells - very few features expressed, likely another class of dying cells
# but either way not super informative and tends to screw up trajectory analyses later on

sce$nexprs <- nexprs(sce) 

discard.nexpres.madlower <- isOutlier(sce$nexprs, type="lower", nmads=1)

sce <- sce[,!discard.nexpres.madlower]

metadata(sce)$n_too_few_genes_exprs <- sum(discard.nexpres.madlower)


# ------------------------------------------------------------------------------
# renormalize (previously used simple norm for scbdlfinder), deconvolution norm for actual analysis

sce$sizeFactor <- NULL
logcounts(sce) <- NULL

clust.pre <- quickCluster(sce)

# occasionally get warning about negative size factors - use positive=T to enforce positivity
sce <- computeSumFactors(sce,cluster=clust.pre, min.mean=0.2, positive=T)

sce <- logNormCounts(sce)

# ------------------------------------------------------------------------------
# variable feature selection

dec <- modelGeneVar(sce)
hvg.var <- getTopHVGs(dec, prop = 0.1)

tes <- rowData(sce)[rowData(sce)$gene_biotype %in% "repeat_element",] |> rownames()

hvtes <- hvg.var[hvg.var %in% tes]

#& 
chosen <- hvg.var[!hvg.var %in% tes & !hvg.var %in% subset(rowData(sce),grepl("^mt-",gene_name))$gene_name]

rowSubset(sce) <- chosen
metadata(sce)$highly.variable.tes <- hvtes
metadata(sce)$highly.variable.genes <- chosen

# ------------------------------------------------------------------------------
# PCA - without influence of mito genes or TEs
sce <- denoisePCA(sce, dec, subset.row=chosen)

g_pc_elbow <- tibble(percent.var = attr(reducedDim(sce), "percentVar")) |>
  mutate(PC = row_number()) |>
  ggplot(aes(PC, percent.var)) +
  geom_point() +
  geom_vline(xintercept = ncol(reducedDim(sce,"PCA")), color="red")


# ------------------------------------------------------------------------------
# UMAP

set.seed(2)
sce <- runUMAP(sce, dimred="PCA")
set.seed(2)
sce <- runTSNE(sce, dimred="PCA")

#plotReducedDim(sce, dimred = "UMAP", colour_by = "subsets_MT_percent",text_by = "x")

# ------------------------------------------------------------------------------
# coarse clustering
set.seed(20)
clust.sweep <- clusterSweep(reducedDim(sce, "PCA"), 
                            NNGraphParam(), 
                            k=as.integer(c(15, 25, 35, 45, 55, 65, 75, 85)),
                            cluster.fun=c("louvain", "walktrap", "infomap","leiden"),
                            BPPARAM=BiocParallel::MulticoreParam(8))

sweep.df <- as_tibble(clust.sweep$parameters,rownames="paramset")

sweep.df <- full_join(sweep.df, enframe(as.list(clust.sweep$clusters),name = "paramset", value = "labels"), by="paramset", relationship="one-to-one")

sweep.df <- sweep.df |>
  mutate(nclust = map_int(labels, ~length(unique(.x))))

sweep.df <- sweep.df |>
  mutate(mean.sil.width= map_dbl(labels,~mean(approxSilhouette(reducedDim(sce,"PCA"),.x)$width)))

sweep.df <- sweep.df |>
  mutate(wcss= map_dbl(labels,~sum(clusterRMSD(reducedDim(sce), .x, sum=TRUE), na.rm=TRUE)))

g_coarse_clustering_sweep <- sweep.df |>
  dplyr::select(-labels) |>
  pivot_longer(-c(paramset,k,cluster.fun),names_to = "metric", values_to = "value") |>
  ggplot(aes(k,value,color=cluster.fun)) +
  facet_wrap(~metric,scales = "free") +
  geom_line()


set.seed(20)
nn.clust <- clusterCells(sce, use.dimred="PCA", full=TRUE,
                         BLUSPARAM=NNGraphParam(cluster.fun="leiden",k =35))

colLabels(sce) <- nn.clust$clusters

#plotReducedDim(sce, dimred = "UMAP", colour_by = "label",text_by = "label") + paletteer::scale_color_paletteer_d("ggsci::default_igv")
#plotReducedDim(sce, dimred = "UMAP", colour_by = "x",text_by = "x") + paletteer::scale_color_paletteer_d("ggsci::default_igv")

sil.approx <- approxSilhouette(reducedDim(sce, "PCA"), clusters=colLabels(sce)) |>
  as_tibble(rownames = "Barcode")

sil.approx$closest <- factor(ifelse(sil.approx$width > 0, colLabels(sce), sil.approx$other))
sil.approx$cluster <- colLabels(sce)

g_silhouette <- ggplot(sil.approx, aes(x=cluster, y=width, colour=closest)) +
  geom_jitter()

# ------------------------------------------------------------------------------
# subclustering
# this takes quite a while to run.
set.seed(2)
subcluster.out <- quickSubCluster(sce, groups=sce$label,simplify=T,
                                  prepFUN=function(x) { # Preparing the subsetted SCE for clustering.
                                    dec <- modelGeneVar(x)
                                    input <- denoisePCA(x, technical=dec,
                                                        subset.row=getTopHVGs(dec, prop=0.1),
                                                        BSPARAM=BiocSingular::IrlbaParam())
                                  },
                                  clusterFUN=function(x) { # Performing the subclustering in the subset.
                                    g <- buildSNNGraph(x, use.dimred="PCA", k=55)
                                    igraph::cluster_leiden(g)$membership
                                  }
)

sce$sublabel <- subcluster.out


# ------------------------------------------------------------------------------
# cluster naming by vote - cells used by umich collaborators get a vote
colData(sce) <- colData(sce) |>
  as_tibble(rownames="cell") |>
  dplyr::count(sublabel, x) |>
  filter(!is.na(x)) |>
  group_by(sublabel) |>
  slice_max(n=1, n) |>
  ungroup() |> 
  mutate(celltype = str_extract(x,"Undifferentiated-SSC|Differentiated-SSC|Preleptotene-SPC|Myoid|TCF21\\+|sertoli|Leydig|Macrophage|Endothelial|11")) |>
  dplyr::select(sublabel, celltype) |>
  left_join(x=as_tibble(colData(sce),rownames="cell"),y=_, by="sublabel") |>
  column_to_rownames("cell") |>
  DataFrame()

sce$sublabel2 <- paste(sce$sublabel,sce$celltype,sep="/")

colData(sce) <- colData(sce) |>
  as_tibble(rownames="cell") |>
  dplyr::count(label, x) |>
  filter(!is.na(x)) |>
  group_by(label) |>
  slice_max(n=1, n) |>
  ungroup() |> 
  mutate(macro_celltype = str_extract(x,"SSC|SPC|Myoid|TCF21\\+|sertoli|Leydig|Macrophage|Endothelial|11")) |>
  dplyr::select(label, macro_celltype) |>
  left_join(x=as_tibble(colData(sce),rownames="cell"),y=_, by="label") |>
  column_to_rownames("cell") |>
  DataFrame()

sce$label2 <- paste(sce$label,sce$macro_celltype,sep="/")

plotReducedDim(sce, dimred="TSNE",colour_by = "celltype")

# ------------------------------------------------------------------------------
# subset by variously useful grouping
#https://www.nature.com/articles/s41467-021-24130-8#:~:text=In%20vivo%2C%20fetal%20TCF21lin,well%20as%20in%20normal%20aging.
macro_grps <- list(sertoli = "sertoli",
                   germ_cell = c("SSC","SPC"),
                   somatic_non_macro = c("TCF21+","Leydig","sertoli","Myoid"),
                   macrophage = "Macrophage" # may encompass endothelial, depending on clustering
                   )

sce_subsets <- map(macro_grps, ~{sce[,sce$macro_celltype %in% .x]})
  
# ------------------------------------------------------------------------------
# export
write_rds(sce, snakemake@output$rds)
write_rds(sce_mito_warning, snakemake@output$sce_mito_warning)
write_rds(g_dcx, snakemake@output$g_dcx)
write_rds(g_mito_cuts, snakemake@output$g_mito_cuts)
write_rds(g_pc_elbow, snakemake@output$g_pc_elbow)
write_rds(g_coarse_clustering_sweep, snakemake@output$g_coarse_clustering_sweep)
write_rds(nn.clust, snakemake@output$nn_clust)
write_rds(g_silhouette, snakemake@output$g_silhouette)


write_rds(sce_subsets$sertoli, snakemake@output$sce_sertoli)
write_rds(sce_subsets$germ_cell, snakemake@output$sce_germ_cell)
write_rds(sce_subsets$somatic_non_macro, snakemake@output$sce_somatic_non_macro)
write_rds(sce_subsets$macrophage, snakemake@output$sce_macrophage)
