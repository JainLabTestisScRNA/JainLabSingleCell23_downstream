library(tidyverse)
library(scater)
library(scran)
library(scuttle)

theme_set(theme_classic())

de_by <-  ifelse(exists("snakemake"),
                 snakemake@params$de_by,
                 "label")

# regex for column used for differential expression (currently 'label' for GC clusters,"celltype" for overall)
de_by <- sprintf("^%s$",de_by)

sce_fl <- ifelse(exists("snakemake"),
                 snakemake@input$sce,
                 "results/germ_cells/adult.sce.germ_cell.both_genotypes.subclustered.reintegrated.rds")

sce <- read_rds(sce_fl)

# ------------------------------------------------------------------------------
# get de results, filtering and ordering by celltypes we want to examine
de_fl <- ifelse(exists("snakemake"),
                snakemake@input$de,
                "results/differential_expression/adult.tbl.germ_cell.reprocessed.de.tsv.gz")

de <- read_tsv(de_fl)

# handle ordering whther using broad anno or numbered germ clusters
if ("Spermatogonia" %in% de$celltype) {
  celltype_ord <- c("Spermatogonia","Spermatocyte","RoundSpermatid","Elongating") 
} else {
  celltype_ord <- unique(de$celltype)
  celltype_ord <- celltype_ord[str_extract(celltype_ord,"\\d+") |> as.integer() |> order()]
}

de <- de |>
  filter(celltype %in% celltype_ord) |>
  mutate(celltype = fct_relevel(celltype,rev(celltype_ord)))

# ------------------------------------------------------------------------------
# get repeat type classiciations and create a df for annotating TE features
classification_fl <- 
  ifelse(exists("snakemake"),
         snakemake@input$classification,
         "data/dfam_classification.parsed.txt")

classification <- read_tsv(classification_fl) |>
  dplyr::select(dfam_name,classification)

# ------------------------------------------------------------------------------
# heatmap of fold changes

fc <- filter(de,!str_detect(feature,"ENSMUSG")) |>
  group_by(feature) |>
  filter(any(FDR<=0.05)) |>
  ungroup() |>
  mutate(logFC = if_else(FDR > 0.05,NA,logFC)) |>
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
  scale_fill_gradient2(breaks=scales::pretty_breaks(n=4)(-max(abs(fc$logFC),na.rm = T):max(abs(fc$logFC),na.rm = T)),limits=c(-max(abs(fc$logFC),na.rm = T),max(abs(fc$logFC),na.rm = T)),low = "blue",high="red") +
  theme(axis.text.x = element_text(angle=45,hjust=1))

ggsave(snakemake@output$pdf_foldchange,g1,height = 4,width = snakemake@params$width)


g2 <- fc |>
  #filter(celltype %in% coi) |>
  mutate(classification = factor(classification)) |>
  mutate(feature = fct_reorder(feature,logFC)) |>
  ggplot(aes(feature,celltype,fill=-log10(FDR))) +
  geom_tile() +
  facet_grid(.~classification,scales="free",space = "free") +
  scale_fill_gradient(low = "white",high="darkgreen") +
  theme(axis.text.x = element_text(angle=45,hjust=1))

ggsave(snakemake@output$pdf_fdr,g2,height = 4,width = snakemake@params$width)

write_tsv(fc,snakemake@output$tsv)