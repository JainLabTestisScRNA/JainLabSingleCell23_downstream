library(scater)
library(scran)
library(scuttle)
library(tidyverse)
library(corrr)

theme_set(theme_classic())

fl_mut <- ifelse(exists("snakemake"),snakemake@input$mut,"results/germ_cells/adult.sce.integrated.clustered.celltypes.germ_cell.mut.reprocessed.rds")
sce_mut <- read_rds(fl_mut)

fl_wt <- ifelse(exists("snakemake"),snakemake@input$wt,"results/germ_cells/adult.sce.integrated.clustered.celltypes.germ_cell.wt.reprocessed.rds")
sce_wt <- read_rds(fl_wt)


#goi <- readxl::read_xlsx("data/green2018_supp3.xlsx",skip = 5) |> pull(gene)

# get correlations green et al clusters and our clusters
gc12_mat <- readxl::read_xlsx("data/GSE112393_MergedAdultMouseST25_12GermCellClusters_AllGeneExp.xlsx",skip = 4) |>
  dplyr::rename(gene_name ="Cluster") |>
  filter(!row_number() %in% c(1,2)) |>
  inner_join(as_tibble(rowData(sce_wt)[,c("gene_id","gene_name")])) |>
  #filter(gene_name %in% goi) |>
  dplyr::select(-gene_name) |>
  dplyr::relocate(gene_id)


# ------------------------------------------------------------------------------
# create pseudobulk for labels, not considering genotype
pse_mut <- aggregateAcrossCells(sce_mut,sce_mut$label)
pse_mut <- computeSumFactors(pse_mut)
pse_mut <- logNormCounts(pse_mut)
colnames(pse_mut) <- paste0(colnames(pse_mut),"_MUT")

pse_wt <- aggregateAcrossCells(sce_wt,sce_wt$label)
pse_wt <- computeSumFactors(pse_wt)
pse_wt <- logNormCounts(pse_wt)
colnames(pse_wt) <- paste0(colnames(pse_wt),"_WT")

c_df <- as_tibble(logcounts(pse_wt),rownames="gene_id") |>
  inner_join(as_tibble(logcounts(pse_mut),rownames="gene_id"),by="gene_id") |>
  inner_join(gc12_mat,by="gene_id")  |>
  corrr::correlate(diagonal=1)

write_tsv(c_df,snakemake@output$mat)


# plot wt vs green et al
pdf(snakemake@output$wt_vs_green2018)
c_df |>
  filter(str_detect(term,"WT")) |>
  dplyr::select(c(term,contains("GC"))) |>
  mutate(term = fct_reorder(term,as.integer(str_extract(term,"\\d+")))) |>
  arrange(term) |>
  column_to_rownames("term") |>
  pheatmap::pheatmap(cluster_cols = F,cluster_rows = F)
dev.off()

# plot mut vs green et al
pdf(snakemake@output$mut_vs_green2018)
c_df |>
  filter(str_detect(term,"MUT")) |>
  dplyr::select(c(term,contains("GC"))) |>
  mutate(term = fct_reorder(term,as.integer(str_extract(term,"\\d+")))) |>
  arrange(term) |>
  column_to_rownames("term") |>
  pheatmap::pheatmap(cluster_cols = F,cluster_rows = F)
dev.off()


# plot separately processed mut vs wt, not considering genotype
pdf(snakemake@output$mut_vs_wt)
x <- c_df |>
  filter(str_detect(term,"MUT")) |>
  dplyr::select(c(term,contains("WT"))) |>
  mutate(term = fct_reorder(term,as.integer(str_extract(term,"\\d+")))) |>
  arrange(term) |>
  column_to_rownames("term")

x <- x[,order(as.integer(str_extract(colnames(x),"\\d+")))]

pheatmap::pheatmap(x, cluster_cols = F,cluster_rows = F)
dev.off()
