library(scran)
library(scuttle)
library(scater)
library(tidyverse)
library(patchwork)
library(slingshot)

de <- read_tsv("results/differential_expression/adult.tbl.pseudobulk.de.tsv.gz")

filter(de,!str_detect(feature,"ENSMUSG") & FDR < 0.05) |>
  ggplot(aes(celltype,logFC)) +
  geom_boxplot() +
  geom_jitter()






scemw <- read_rds("results/single-cell-preproc/preprocessed/adult_9646_combined.cellranger.sce.mito_warning.rds")
sce <- read_rds("results/single-cell-preproc/preprocessed/adult_9646_combined.cellranger.sce.passing_cells.rds")

any(colnames(scemw) %in% colnames(sce))

df <- makePerCellDF(sce) |>
  as_tibble() |>
  dplyr::select(contains("subsets"),mito_warning)

dfmw <-  makePerCellDF(scemw) |>
  as_tibble() |>
  dplyr::select(contains("subsets"),mito_warning)


bind_rows(df,dfmw) |>
  ggplot(aes(subsets_MT_percent,subsets_TE_percent,color=mito_warning)) +
  geom_point() +
  scale_y_log10()

germ_sce <- read_rds("results/trajectory/adult.cellranger.sce.germ_cell.trajectory.rds")
germ_sce_exportable <- germ_sce
germ_sce_exportable$slingshot <- NULL

embedded <- embedCurves(germ_sce, "UMAP")
embedded <- slingCurves(embedded)[[1]] # only 1 path.
embedded <- data.frame(embedded$s[embedded$ord,])

get_ps_umap <- \(x) plotUMAP(x, colour_by="celltype") +
  geom_path(data=embedded[seq(nrow(embedded)-50),], aes(x=Dim.1, y=Dim.2), size=1.2, arrow = arrow(length=unit(0.30,"cm"), ends="last", type = "closed")) +
  labs(title=unique(x$genotype))

get_ps_umap(germ_sce[,germ_sce$genotype==c("WT")]) +
get_ps_umap(germ_sce[,germ_sce$genotype==c("MUT")])


#https://github.com/kstreet13/slingshot/issues/42
makePerCellDF(germ_sce_exportable,features = c("L1MdA_I_orf2","L1MdA_I_5end","L1MdA_I_3end")) |>
  as_tibble(rownames="cell") |>
  arrange(L1MdA_I_orf2) |>
  ggplot(aes(slingPseudotime_1,fct_reorder(label2,slingPseudotime_1,.na_rm = T,median),color=L1MdA_I_orf2)) +
  geom_jitter(width = 0.1) +
  facet_wrap(~batch + genotype)

# ------ propoprtion - see greenbaum 2006 for inspiration for this plot
sce <- read_rds("results/single-cell-preproc/integrated/adult.cellranger.sce.integrated.rds")

df <- makePerCellDF(sce) |>
  as_tibble(rownames="cell")


#https://stackoverflow.com/questions/43769234/r-grouped-barplot-of-proportions-of-character-variables-with-error-bars
df |>
  group_by(batch,celltype, genotype) |>
  tally() |>
  group_by(batch,genotype) |>
  mutate(n_tot=sum(n)) |>
  ungroup() |>
  mutate(prop=n/n_tot,
         error = sqrt((prop*(1-prop))/n)) |>
  filter(celltype %in% c("Spermatogonia","Spermatocyte","Elongating","RoundSpermatid")) |>
  mutate(celltype=fct_relevel(celltype, c("Spermatogonia","Spermatocyte","Elongating","RoundSpermatid"))) |>
  mutate(ymin=prop-error) |>
  mutate(ymin=map_dbl(ymin,max,0)) |>
  ggplot(aes(celltype,prop,fill=genotype)) +
  facet_wrap(~batch) +
  geom_errorbar(aes(ymin = ymin,ymax = prop + error),width=0.75,
                position = position_dodge(0.9)) +
  geom_col(position = "dodge") +
  ylab("proportion post-filtering cells") +
  theme_classic() +
  scale_fill_grey() +
  theme(axis.text.x = element_text(angle=45, hjust=1))


# --------------- dysreg in adjult vs juvenile

adult_dysreg <- names(Biostrings::readDNAStringSet("results/tradeseq/adult.cellranger.germ_cell.dysregulated_tes.fasta"))

juv_dysreg <- names(Biostrings::readDNAStringSet("results/tradeseq/juvenile_13d_wt_null.cellranger.germ_cell.dysregulated_tes.fasta"))

dysreg <- list(adult=adult_dysreg, juvenile=juv_dysreg) |>
  enframe() |>
  unnest(value) |>
  mutate(dysregulated=T) |>
  pivot_wider(names_from = name, values_from = dysregulated,values_fill = F)

filter(dysreg,str_detect(value,"^L1")) |>
  gt::gt()

filter(dysreg,!str_detect(value,"^L1")) |>
  gt::gt()

intersect(adult_dysreg, juv_dysreg)


# --- plot dazl and tex19.1 expression

germ_sce |>
  plotExpression(features = c("Dazl","Tex14","Tex19.1","Phf7"), swap_rownames = "gene_name",colour_by = "genotype", other_fields = "celltype") +
  facet_wrap(~ celltype+colour_by,ncol=2)

germ_sce |>
  plotExpression(features = c("Eif4e","Phf7"), swap_rownames = "gene_name",colour_by = "genotype", other_fields = "celltype") +
  facet_wrap(~ celltype+colour_by,ncol=2)

plotReducedDim(germ_sce,dimred="UMAP",colour_by = c("Tex14"),
               other_fields = c("genotype"),
               swap_rownames = "gene_name") +
  facet_wrap(~genotype)

plotReducedDim(germ_sce,dimred="UMAP",colour_by = c("Dazl"),
               other_fields = c("genotype"),
               swap_rownames = "gene_name") +
  facet_wrap(~genotype)

plotReducedDim(germ_sce,dimred="UMAP",colour_by = c("Tex19.1"),
               other_fields = c("genotype"),
               swap_rownames = "gene_name") +
  facet_wrap(~genotype)

plotReducedDim(germ_sce,dimred="UMAP",colour_by = c("Ddx4"),
               other_fields = c("genotype"),
               swap_rownames = "gene_name") +
  facet_wrap(~genotype)

plotReducedDim(germ_sce,dimred="UMAP",colour_by = c("Dnmt3c"),
               other_fields = c("genotype"),
               swap_rownames = "gene_name") +
  facet_wrap(~genotype)

plotReducedDim(germ_sce,dimred="UMAP",colour_by = c("Phf7"),
               other_fields = c("genotype"),
               swap_rownames = "gene_name") +
  facet_wrap(~genotype)


deg <- read_rds("results/tradeseq/adult.cellranger.germ_cell.deg.rds")

deg |>
  filter(padj < 0.0001) |>
  arrange(-abs(waldStat)) |>
  head(100) |>
  print(n=Inf)

res <- read_rds("results/tradeseq/adult.cellranger.sce.germ_cell.tradeseq.rds")
res$slingshot <- NULL
res
library(tradeSeq)
tradeSeq::plotSmoothers(res,counts=logcounts(res),gene="L1MdA_I_5end")
