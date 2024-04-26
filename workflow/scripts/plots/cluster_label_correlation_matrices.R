Sys.setenv(R_PROFILE=".Rprofile")
source(Sys.getenv("R_PROFILE"))

library(scater)
library(scran)
library(scuttle)
library(tidyverse)
library(corrr)
library(RColorBrewer)
theme_set(theme_classic())

fl <- ifelse(exists("snakemake"),snakemake@input$sce,"results/germ_cells/adult.sce.germ_cell.both_genotypes.subclustered.reintegrated.rds")
sce <- read_rds(fl)


#goi <- readxl::read_xlsx("data/green2018_supp3.xlsx",skip = 5) |> pull(gene)


coul <- brewer.pal(11, "PiYG") 
coul <- colorRampPalette(coul)(100)

# ------------------------------------------------------------------------------
# create pseudobulk for labels, not considering genotype
pse <- aggregateAcrossCells(sce,paste0(sce$label,"_bothGenotypes"))
pse <- computeSumFactors(pse)
pse <- logNormCounts(pse)

# pseudobulk for genotype/label combos
pse.gt <- aggregateAcrossCells(sce,paste(sce$label,sce$genotype,sep="_"))
pse.gt <- computeLibraryFactors(pse.gt)
pse.gt <- logNormCounts(pse.gt)


# get correlations green et al clusters and our clusters
gc12_mat <- readxl::read_xlsx("data/GSE112393_MergedAdultMouseST25_12GermCellClusters_AllGeneExp.xlsx",skip = 4) |>
  dplyr::rename(gene_name ="Cluster") |>
  filter(!row_number() %in% c(1,2)) |>
  inner_join(as_tibble(makePerFeatureDF(pse)[,c("gene_id","gene_name")])) |>
  #filter(gene_name %in% goi) |>
  dplyr::select(-gene_name) |>
  dplyr::relocate(gene_id)

c_df <-  as_tibble(logcounts(pse), rownames="gene_id") |>
  inner_join(as_tibble(logcounts(pse.gt), rownames="gene_id"),by="gene_id") |>
  inner_join(gc12_mat,by="gene_id")  |>
  corrr::correlate(diagonal=1,method = "pearson")

write_tsv(c_df,snakemake@output$mat)


# plot us vs green et al
pdf(snakemake@output$jain_vs_green2018)
c_df |>
  filter(str_detect(term,"both")) |>
  dplyr::select(c(term,contains("GC"))) |>
  mutate(term = fct_reorder(term,as.integer(str_extract(term,"\\d+")))) |>
  arrange(term) |>
  column_to_rownames("term") |>
  pheatmap::pheatmap(cluster_cols = F,cluster_rows = F,color = coul)
dev.off()

# plot green et al vs green et al

pdf(snakemake@output$green2018_vs_green2018)
c_df |>
  filter(str_detect(term,"GC")) |>
  dplyr::select(c(term,contains("GC"))) |>
  mutate(term = fct_reorder(term,as.integer(str_extract(term,"\\d+")))) |>
  arrange(term) |>
  column_to_rownames("term") |>
  pheatmap::pheatmap(cluster_cols = F,cluster_rows = F,color = coul)
dev.off()


# plot us vs us, not considering genotype
pdf(snakemake@output$jain_vs_jain)
x <- c_df |>
  filter(str_detect(term,"both")) |>
  dplyr::select(c(term,contains("both"))) |>
  mutate(term = fct_reorder(term,as.integer(str_extract(term,"\\d+")))) |>
  arrange(term) |>
  column_to_rownames("term")
x <- x[,rownames(x)]

pheatmap::pheatmap(x,cluster_cols = F,cluster_rows = F,color = coul)
dev.off()

# ------------------------------------------------------------------------------
# recompute correlations for only our data, considering genotype
# ------------------------------------------------------------------------------

# plot us vs green et al
pdf(snakemake@output$mut_vs_green2018)
c_df |>
  filter(str_detect(term,"MUT")) |>
  dplyr::select(c(term,contains("GC"))) |>
  mutate(term = fct_reorder(term,as.integer(str_extract(term,"\\d+")))) |>
  arrange(term) |>
  column_to_rownames("term") |>
  pheatmap::pheatmap(cluster_cols = F,cluster_rows = F,color = coul)
dev.off()

# plot us vs green et al
pdf(snakemake@output$wt_vs_green2018)
c_df |>
  filter(str_detect(term,"WT")) |>
  dplyr::select(c(term,contains("GC"))) |>
  mutate(term = fct_reorder(term,as.integer(str_extract(term,"\\d+")))) |>
  arrange(term) |>
  column_to_rownames("term") |>
  pheatmap::pheatmap(cluster_cols = F,cluster_rows = F,color = coul)
dev.off()


pdf(snakemake@output$mut_vs_wt)
x <- c_df |>
  filter(str_detect(term,"MUT")) |>
  dplyr::select(c(term,contains("WT"))) |>
  dplyr::select(term,matches("_WT")) |>
  mutate(term = fct_reorder(term,as.integer(str_extract(term,"\\d+")))) |>
  arrange(term) |>
  column_to_rownames("term")

pheatmap::pheatmap(x[,order(as.integer(str_extract(colnames(x),"\\d+")))],cluster_cols = F,cluster_rows = F,color = coul)
dev.off()

pdf(snakemake@output$wt_vs_wt)
x <- c_df |>
  filter(str_detect(term,"WT")) |>
  dplyr::select(c(term,contains("WT"))) |>
  mutate(term = fct_reorder(term,as.integer(str_extract(term,"\\d+")))) |>
  arrange(term) |>
  column_to_rownames("term")
pheatmap::pheatmap(x[,order(as.integer(str_extract(colnames(x),"\\d+")))],cluster_cols = F,cluster_rows = F,color = coul)
dev.off()

pdf(snakemake@output$mut_vs_mut)
x <- c_df |>
  filter(str_detect(term,"MUT")) |>
  #filter(!term %in% c("cl10_MUT","cl11_MUT","cl12_MUT")) |>
  dplyr::select(c(term,contains("MUT"))) |>
  mutate(term = fct_reorder(term,as.integer(str_extract(term,"\\d+")))) |>
  arrange(term) |>
  column_to_rownames("term")
pheatmap::pheatmap(x[,order(as.integer(str_extract(colnames(x),"\\d+")))],cluster_cols = F,cluster_rows = F,color = coul)
dev.off()
