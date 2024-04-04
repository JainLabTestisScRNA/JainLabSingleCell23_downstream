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

de_fl <- ifelse(exists("snakemake"),
                snakemake@input$de,
                "results/differential_expression/adult.tbl.pseudobulk.de.tsv.gz")

de <- read_tsv(de_fl)

classification_fl <- 
  ifelse(exists("snakemake"),
         snakemake@input$classification,
         "data/dfam_classification.parsed.txt")

classification <- read_tsv(classification_fl) |>
  dplyr::select(dfam_name,classification)


df <- filter(de,!str_detect(feature,"ENSMUSG")) |>
  left_join(classification,by=c(feature="dfam_name")) |>
  filter(FDR < 0.05) |>
  #filter(classification == "LINE") |>
  dplyr::select(feature,classification) |>
  arrange(classification) |>
  distinct()

# ------------------------------------------------------------------------------
# heatmap of fold changes
celltype_ord <- c("Spermatogonia","Spermatocyte","RoundSpermatid","Elongating",
                  "Sertoli","InnateLymphoid","Macrophage","Endothelial") 

fc <- filter(de,!str_detect(feature,"ENSMUSG")) |>
  group_by(feature) |>
  filter(any(FDR<0.05)) |>
  ungroup() |>
  left_join(classification,by=c(feature="dfam_name")) |>
  #filter(classification == "LINE") |>
  dplyr::select(feature,classification,celltype,logFC) |>
  arrange(classification) |>
  distinct() |>
  mutate(celltype = fct_relevel(celltype,rev(celltype_ord)))

g1 <- fc |>
  #filter(celltype %in% coi) |>
  mutate(classification = factor(classification)) |>
  mutate(feature = fct_reorder(feature,logFC)) |>
  ggplot(aes(feature,celltype,fill=logFC)) +
  geom_tile() +
  facet_grid(.~classification,scales="free",space = "free") +
  scale_fill_gradient2(breaks=scales::pretty_breaks(n=4)(-1.5:1.5),limits=c(-1.5,1.5),low = "blue",high="red") +
  theme(axis.text.x = element_text(angle=45,hjust=1))

ggsave(snakemake@output$pdf_foldchange,g1,height = 4,width = 8)

# ------------------------------------------------------------------------------
# heatmap of per cell expression

pdf(snakemake@output$pdf)
coi <- c("Spermatogonia","Spermatocyte","RoundSpermatid","Elongating")
sce.g <- sce[df$feature,sce$celltype %in% coi]
sce.g$plotgroup <- paste(sce.g$celltype,sce.g$genotype)

rowData(sce.g)$classification <- rowData(sce.g)$gene_symbol |> map_chr(~deframe(df)[[.x]])

plotHeatmap(sce.g,scale = T,center=T,zlim = c(-1.5,1.5),
            features = df$feature,
            color_rows_by = "classification",
            order_rows_by = "classification",
            order_columns_by = c("genotype","celltype","batch","nexprs"),
            color_columns_by = c("genotype","batch","celltype"))


# ------------------------------------------------------------------------------
# heatmap of averaged per group expression
plotGroupedHeatmap(sce.g,
                   scale = T,center=T,zlim = c(-1.5,1.5),
                   color_rows_by = "classification",
                   block = c("batch"),group = c("plotgroup"),features = df$feature)

dev.off()

write_tsv(fc,snakemake@output$tsv)