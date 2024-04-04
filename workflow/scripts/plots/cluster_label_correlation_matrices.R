library(scater)
library(scran)
library(scuttle)
library(tidyverse)
library(corrr)

theme_set(theme_classic())

fl <- ifelse(exists("snakemake"),snakemake@input$sce,"results/germ_cells/adult.sce.integrated.clustered.celltypes.germ_cell.reprocessed.rds")
sce <- read_rds(fl)


goi <- readxl::read_xlsx("data/green2018_supp3.xlsx",skip = 5) |> pull(gene)

# ------------------------------------------------------------------------------
# create pseudobulk for labels, not considering genotype
pse <- aggregateAcrossCells(sce,sce$label)
pse <- computeSumFactors(pse)
pse <- logNormCounts(pse)

# get correlations green et al clusters and our clusters
gc12_mat <- readxl::read_xlsx("data/GSE112393_MergedAdultMouseST25_12GermCellClusters_AllGeneExp.xlsx",skip = 4) |>
  dplyr::rename(gene_name ="Cluster") |>
  filter(!row_number() %in% c(1,2)) |>
  inner_join(as_tibble(makePerFeatureDF(pse)[,c("gene_id","gene_name")])) |>
  filter(gene_name %in% goi) |>
  dplyr::select(-gene_name) |>
  dplyr::relocate(gene_id)

c_df <- logcounts(pse) |> as_tibble(rownames="gene_id") |>
  inner_join(gc12_mat,by="gene_id")  |>
  corrr::correlate(diagonal=1)

write_tsv(c_df,snakemake@output$no_genotype)


# plot us vs green et al
pdf(snakemake@output$jain_vs_green2018)
c_df |>
  filter(str_detect(term,"cl")) |>
  dplyr::select(c(term,contains("GC"))) |>
  mutate(term = fct_reorder(term,as.integer(str_extract(term,"\\d+")))) |>
  arrange(term) |>
  column_to_rownames("term") |>
  pheatmap::pheatmap(cluster_cols = F,cluster_rows = F)
dev.off()

# plot green et al vs green et al

pdf(snakemake@output$green2018_vs_green2018)
c_df |>
  filter(str_detect(term,"GC")) |>
  dplyr::select(c(term,contains("GC"))) |>
  mutate(term = fct_reorder(term,as.integer(str_extract(term,"\\d+")))) |>
  arrange(term) |>
  column_to_rownames("term") |>
  pheatmap::pheatmap(cluster_cols = F,cluster_rows = F)
dev.off()


# plot us vs us, not considering genotype
pdf(snakemake@output$jain_vs_jain)
c_df |>
  filter(str_detect(term,"cl")) |>
  dplyr::select(c(term,contains("cl"))) |>
  dplyr::select(term,cl1,cl2,cl3,cl4,cl5,cl6,cl7,cl8,cl9,cl10,cl11,cl12) |>
  mutate(term = fct_reorder(term,as.integer(str_extract(term,"\\d+")))) |>
  arrange(term) |>
  column_to_rownames("term") |>
  pheatmap::pheatmap(cluster_cols = F,cluster_rows = F)
dev.off()

# ------------------------------------------------------------------------------
# recompute correlations for only our data, considering genotype
# ------------------------------------------------------------------------------

# pseudobulk for genotype/label combos
pse <- aggregateAcrossCells(sce,paste(sce$label,sce$genotype,sep="_"))
pse <- computeLibraryFactors(pse)
pse <- logNormCounts(pse)

c_df <- logcounts(pse) |> as_tibble(rownames="gene_id") |>
  corrr::correlate(diagonal=1)

write_tsv(c_df,snakemake@output$genotype)

# standardize the colors for these plots
palette_len <- 10
colors <- colorRampPalette( c("blue","white","red"))(palette_len)
#colors <- viridis::viridis(palette_len,direction = 1,begin = 1,end=0.1)

breaks <- seq(min(c_df[,-1]), max(c_df[-1]), length.out =palette_len)

pdf(snakemake@output$mut_vs_wt)
c_df |>
  filter(str_detect(term,"MUT")) |>
  dplyr::select(c(term,contains("WT"))) |>
  dplyr::select(term,cl1_WT,cl2_WT,cl3_WT,cl4_WT,cl5_WT,cl6_WT,cl7_WT,cl8_WT,cl9_WT,cl10_WT,cl11_WT,cl12_WT) |>
  mutate(term = fct_reorder(term,as.integer(str_extract(term,"\\d+")))) |>
  arrange(term) |>
  column_to_rownames("term") |>
  pheatmap::pheatmap(cluster_cols = F,cluster_rows = F,breaks = breaks, color = colors)
dev.off()

pdf(snakemake@output$wt_vs_wt)
c_df |>
  filter(str_detect(term,"WT")) |>
  dplyr::select(c(term,contains("WT"))) |>
  dplyr::select(term,cl1_WT,cl2_WT,cl3_WT,cl4_WT,cl5_WT,cl6_WT,cl7_WT,cl8_WT,cl9_WT,cl10_WT,cl11_WT,cl12_WT) |>
  mutate(term = fct_reorder(term,as.integer(str_extract(term,"\\d+")))) |>
  arrange(term) |>
  column_to_rownames("term") |>
  pheatmap::pheatmap(cluster_cols = F,cluster_rows = F,breaks = breaks, color = colors)
dev.off()

pdf(snakemake@output$mut_vs_mut)
c_df |>
  filter(str_detect(term,"MUT")) |>
  #filter(!term %in% c("cl10_MUT","cl11_MUT","cl12_MUT")) |>
  dplyr::select(c(term,contains("MUT"))) |>
  dplyr::select(term,cl1_MUT,cl2_MUT,cl3_MUT,cl4_MUT,cl5_MUT,cl6_MUT,cl7_MUT,cl8_MUT,cl9_MUT,cl10_MUT,cl11_MUT,cl12_MUT) |>
  mutate(term = fct_reorder(term,as.integer(str_extract(term,"\\d+")))) |>
  arrange(term) |>
  column_to_rownames("term") |>
  pheatmap::pheatmap(cluster_cols = F,cluster_rows = F,breaks = breaks, color = colors)
dev.off()
