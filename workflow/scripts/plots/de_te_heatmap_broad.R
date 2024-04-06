library(tidyverse)
library(scater)
library(scran)
library(scuttle)

theme_set(theme_classic())

de_by <-  ifelse(exists("snakemake"),
                 snakemake@params$de_by,
                 "celltype")

# regex for column used for differential expression (currently 'label' for GC clusters,"celltype" for overall)
de_by <- sprintf("^%s$",de_by)

sce_fl <- ifelse(exists("snakemake"),
                 snakemake@input$sce,
                 "results/integration/adult.sce.integrated.clustered.celltypes.rds")

sce <- read_rds(sce_fl)

# ------------------------------------------------------------------------------
# get de results, filtering and ordering by celltypes we want to examine
de_fl <- ifelse(exists("snakemake"),
                snakemake@input$de,
                "results/differential_expression/adult.tbl.pseudobulk.de.tsv.gz")

de <- read_tsv(de_fl)

celltype_ord <- c("Spermatogonia","Spermatocyte","RoundSpermatid","Elongating") 
de <- de |>
  filter(celltype %in% celltype_ord) |>
  mutate(celltype = fct_relevel(celltype,rev(celltype_ord)))


# ------------------------------------------------------------------------------
# get repeat type classiciations and creat a df for annotating TE features
classification_fl <- 
  ifelse(exists("snakemake"),
         snakemake@input$classification,
         "data/dfam_classification.parsed.txt")

classification <- read_tsv(classification_fl) |>
  dplyr::select(dfam_name,classification)


#classification <- filter(de,!str_detect(feature,"ENSMUSG")) |>
#  left_join(classification,by=c(feature="dfam_name")) |>
#  filter(FDR < 0.05) |>
#  #filter(classification == "LINE") |>
#  dplyr::select(feature,classification) |>
#  arrange(classification) |>
#  distinct()

# ------------------------------------------------------------------------------
# heatmap of fold changes

fc <- filter(de,!str_detect(feature,"ENSMUSG")) |>
  group_by(feature) |>
  filter(any(FDR<0.05)) |>
  ungroup() |>
  left_join(classification,by=c(feature="dfam_name")) |>
  #filter(classification == "LINE") |>
  dplyr::select(feature,classification,celltype,logFC,FDR) |>
  arrange(classification) |>
  distinct()

g1 <- fc |>
  #filter(celltype %in% coi) |>
  mutate(classification = factor(classification)) |>
  mutate(feature = fct_reorder(feature,logFC)) |>
  ggplot(aes(feature,celltype,fill=logFC)) +
  geom_tile() +
  facet_grid(.~classification,scales="free",space = "free") +
  scale_fill_gradient2(breaks=scales::pretty_breaks(n=4)(-max(abs(fc$logFC)):max(abs(fc$logFC))),limits=c(-max(abs(fc$logFC)),max(abs(fc$logFC))),low = "blue",high="red") +
  theme(axis.text.x = element_text(angle=45,hjust=1))

ggsave(snakemake@output$pdf_foldchange,g1,height = 4,width = 8)


g2 <- fc |>
  #filter(celltype %in% coi) |>
  mutate(classification = factor(classification)) |>
  mutate(feature = fct_reorder(feature,logFC)) |>
  ggplot(aes(feature,celltype,fill=-log10(FDR))) +
  geom_tile() +
  facet_grid(.~classification,scales="free",space = "free") +
  scale_fill_gradient(low = "white",high="darkgreen") +
  theme(axis.text.x = element_text(angle=45,hjust=1))

ggsave(snakemake@output$pdf_fdr,g2,height = 4,width = 8)

write_tsv(fc,snakemake@output$tsv)