library(tidyverse)
library(clusterProfiler)
library(enrichplot)
library(org.Mm.eg.db)
library(scater)

sce_fl <- "results/single-cell-preproc/preprocessed/juvenile_13d_wt_null.cellranger.sce.germ_cell.rds"
sce_fl <- snakemake@input$sce

sce <- read_rds(sce_fl)

# expressed in at least 10 cells
universe <- rownames(sce)[which(nexprs(sce, byrow=T, detection_limit=2) > 10)] |> 
  bitr(fromType = "ENSEMBL", toType = c("ENTREZID"), OrgDb = org.Mm.eg.db) |>
  drop_na() |>
  pull(ENTREZID)

smoothed_fl <- "results/tradeseq/juvenile_13d_wt_null.cellranger.germ_cell.smoothed_cluster_timecourse.rds"
smoothed_fl <- snakemake@input$smoothed
clusters_tsv <- "results/tradeseq/juvenile_13d_wt_null.cellranger.germ_cell.gene_clusters.tsv"
clusters_tsv <- snakemake@input$gene_clusters_tsv

smoothed_df <- read_rds(smoothed_fl)

# get trend of up or down
smoothed_df <- smoothed_df |>
  filter(condition=="mutant") |>
  group_by(feature, gene_name, lineage, cluster) |>
  summarise(lfc = log2(mean(map_dbl(lfcMutant.WT, ~{2^.x}))), .groups = "drop") |> # log2 mean fold change within each cluster
  arrange(lfc)

idf <- smoothed_df |>  
  filter(str_detect(feature,"ENSMUSG")) |>
  pull(feature) |>
  bitr(fromType = "ENSEMBL", toType = c("ENTREZID"), OrgDb = org.Mm.eg.db) |>
  as_tibble() |>
  dplyr::rename(feature="ENSEMBL") |>
  left_join(smoothed_df, y=_, by="feature")

idf <- idf |>
  mutate(direction = if_else(sign(lfc) ==-1,"dn","up"))

forenrich <- filter(idf,str_detect(feature,"ENSMUSG") & !is.na(ENTREZID))

rgo <- clusterProfiler::compareCluster(ENTREZID ~ cluster, data = forenrich, fun="enrichGO", OrgDb= org.Mm.eg.db, ont="ALL", pool=T, universe=universe, readable=T)
rgo.direction <- clusterProfiler::compareCluster(ENTREZID ~ direction*cluster, data = forenrich, fun="enrichGO", OrgDb= org.Mm.eg.db, ont="ALL", pool=T, universe=universe, readable=T)


saveRDS(list(per.cluster =rgo, per.cluster.direction = rgo.direction), snakemake@output$clusterprofiler_rds)

#rgo2 <- setReadable(rgo, 'org.Mm.eg.db', 'ENTREZID')
#rgo.direction2 <- setReadable(rgo.direction, 'org.Mm.eg.db', 'ENTREZID')
#cnetplot(rgo2)
#dotplot(rgo) + theme(axis.text.y = element_text(size=rel(0.7)))
#rgo2.2 <- pairwise_termsim(rgo2)

#p1 <- treeplot(rgo2.2, offset_tiplab = 5, nwords=3)


#emapplot(rgo2.2)
