library(tidyverse)
library(DropletUtils)
library(scater)
library(scran)
library(scuttle)
library(bluster)
library(PCAtools)
library(batchelor)
library(TSCAN)

sce_fls <- Sys.glob("results/germ_cells/adult.sce.germ_cell.wt.*.subclustered.rds")
sce_fls <-  snakemake@input$sces
sces <- map(sce_fls,read_rds)

sces <- map(sces,~{rowSubset(.x) <- NULL;reducedDim(.x) <-  NULL;.x})

all.dec <- lapply(sces, modelGeneVar)

combined.dec <- do.call(combineVar, all.dec)
chosen <- getTopHVGs(combined.dec, prop=0.1)

# don't consider TEs
chosen <- chosen[str_detect(chosen,"ENSMUSG")]

sce <- do.call(cbind,sces)

# ------------------------------------------------------------------------------
# variable feature selection
#dec <- modelGeneVar(sce,block=sce$batch)
#hvg.var <- getTopHVGs(dec, prop = 0.1)
  
#tes <- rowData(sce)[rowData(sce)$gene_biotype %in% "repeat_element" ,] |> rownames()
  
#hvtes <- hvg.var[hvg.var %in% tes]
  
#chosen <- hvg.var[!hvg.var %in% tes & !hvg.var %in% subset(rowData(sce),grepl("^mt-",gene_name))$gene_name]
  
#rowSubset(sce) <- chosen
#metadata(sce)$highly.variable.tes <- hvtes
#metadata(sce)$highly.variable.genes <- chosen
  
# ------------------------------------------------------------------------------
# PCA - without influence of TEs
# ------------------------------------------------------------------------------
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
mnn.out <- fastMNN(sce, k=length(unique(sce$label)),
                     batch = sce$batch, 
                     subset.row=chosen, # no tes involved
                     correct.all = T, 
                     BSPARAM=BiocSingular::RandomParam(deferred = T))
  
reducedDims(sce) <- reducedDims(mnn.out)
assay(sce,"reconstructed") <- assay(mnn.out,"reconstructed")

#plotReducedDim(sce,"corrected",colour_by = "Uchl1",ncomponents = c(1,2),swap_rownames = "gene_name")
#plotReducedDim(sce,"corrected",colour_by = "label",ncomponents = c(1,2),swap_rownames = "gene_name")

pseudo.out <- quickPseudotime(sce,use.dimred="corrected", 
                              use.median = F,
                              dist.method="slingshot",
                              columns=c(1,2),
                              start="1/Spermatogonia",clusters = sce$label,endpoints="1/Elongating")

common.pseudo <- averagePseudotime(pseudo.out$ordering,1)

sce$pseudotime <- common.pseudo

# rename labels by order
label_lookup <- makePerCellDF(sce)[,c("label","celltype","pseudotime")] |>
  as_tibble(rownames = "cell") |>
  mutate(labelNum = str_extract(label,"\\d+")) |>
  group_by(labelNum,label,celltype) |>
  summarise(pseudotime = mean(pseudotime),.groups = "drop") |>
  mutate(newlabel = paste(rank(round(pseudotime,1),ties="min"),celltype,sep="/")) |>
  dplyr::select(label,newlabel) |>
  deframe()

sce$label <- map_chr(sce$label,~{label_lookup[.x]})

plotReducedDim(sce, "corrected",colour_by=I(common.pseudo), 
         text_by="label", text_colour="red") +
  geom_line(data=pseudo.out$connected$corrected, 
            mapping=aes(x=dim1, y=dim2, group=edge))

#plotReducedDim(sce, "corrected",colour_by="label", 
#               text_by="label", text_colour="red") +
#  geom_line(data=pseudo.out$connected$corrected, 
#            mapping=aes(x=dim1, y=dim2, group=edge))

write_rds(sce,snakemake@output$rds)



# ------------------------------------------------------------------------------
# create pseudobulk for labels, not considering genotype
pse <- aggregateAcrossCells(sce,sce$label)
pse <- computeSumFactors(pse)
pse <- logNormCounts(pse)

goi <- readxl::read_xlsx("data/green2018_supp3.xlsx",skip = 5) |> pull(gene)
# get correlations green et al clusters and our clusters
gc12_mat <- readxl::read_xlsx("data/GSE112393_MergedAdultMouseST25_12GermCellClusters_AllGeneExp.xlsx",skip = 4) |>
  dplyr::rename(gene_name ="Cluster") |>
  filter(!row_number() %in% c(1,2)) |>
  inner_join(as_tibble(makePerFeatureDF(pse)[,c("gene_id","gene_name")])) |>
  #filter(gene_name %in% goi) |>
  dplyr::select(-gene_name) |>
  dplyr::relocate(gene_id)

c_df <-  as_tibble(logcounts(pse), rownames="gene_id") |>
  inner_join(gc12_mat,by="gene_id")  |>
  corrr::correlate(diagonal=1,method = "spearman")


# plot us vs green et al
c_df |>
  filter(str_detect(term,"\\/")) |>
  dplyr::select(c(term,contains("GC"))) |>
  mutate(term = fct_reorder(term,as.integer(str_extract(term,"\\d+")))) |>
  arrange(term) |>
  column_to_rownames("term") |>
  pheatmap::pheatmap(cluster_cols = F,cluster_rows = F)



