Sys.setenv(R_PROFILE=".Rprofile")
source(Sys.getenv("R_PROFILE"))

library(tidyverse)
library(scater)
library(scran)
library(scuttle)
library(ComplexHeatmap)
library(RColorBrewer)

de <- read_tsv(ifelse(exists("snakemake"),snakemake@input$de,"results/differential_expression/adult.tbl.germ_cell.reprocessed.de.tsv.gz"))# |>

de_tes <- de |>
  filter(FDR < 0.05 & !str_detect(feature,"ENSMUSG")) |>
  pull(feature) |>
  unique()
#filter(classification == "LINE")

coul <- brewer.pal(11, "PiYG") 
coul <- colorRampPalette(coul)(100)

tes <- read_tsv(ifelse(exists("snakemake"),snakemake@input$classifications,"data/dfam_classification.parsed.txt"))# |>
#filter(classification == "LINE")

te_lkup <- tes |> dplyr::select(dfam_name,classification) |> deframe()

sce <- read_rds(ifelse(exists("snakemake"),
                 snakemake@input$sce,
                 "results/germ_cells/adult.sce.germ_cell.both_genotypes.subclustered.reintegrated.rds"))

sce2 <- aggregateAcrossCells(sce,colData(sce)[,c("label","genotype")],coldata.merge=list(pseudotime=mean),statistics="sum")
sizeFactors(sce2) <- librarySizeFactors(sce2)
sce2 <- logNormCounts(sce2)

# only plot TEs, plus light filtering to avoid uninformative TEs
sce2 <- sce2[rownames(sce2)%in%tes$dfam_name,]
sce2 <- sce2[rowSums(counts(sce2)) > 1,]

# include info about te types
rowData(sce2)$type <- map_chr(rownames(sce2),~{te_lkup[.x]})
rd <- rowData(sce2)[,"type",drop=F] |> as.data.frame()
rd$de <- if_else(rownames(rd) %in% de_tes,"(DE)","")
rd$lab <- sprintf("%s %s",rd$type,rd$de)

assay(sce2,"scaled") <- logcounts(sce2) |> t() |> scale() |> t()

colnames(sce2) <- paste(sce2$label,"_",sce2$genotype)

plot_heat <-  function(feats,rn=F,order_by_cl3=F) {
  # row ord/clustering based on both gt
  rdend <- hclust(dist(assay(sce2,"scaled")[feats,]))
  m <- assay(sce2,"scaled")[feats,]
  if (order_by_cl3) {
    m <- m[order(m[,5],decreasing = T),]
  } else {
    m <- m[rdend$order,]
  }
  Heatmap(m,cluster_rows = F,cluster_columns=F,show_row_names = rn,show_column_names = T,col = rev(coul),name = "Expression\n(row z-scores)",
          row_dend_width = unit(20, "mm"), 
          split=factor(rd[rownames(m),"lab"]),
          border = "black",
          column_labels = sce2$label,column_split = sce2$genotype)
}


pdf(snakemake@output$pdf)

h <- plot_heat(filter(tes,dfam_name %in% rownames(sce2) & classification  %in% c("LINE"))$dfam_name,order_by_cl3 = T)
draw(h,column_title="LINEs")

h <- plot_heat(filter(tes,dfam_name %in% rownames(sce2) & classification  %in% c("LINE"))$dfam_name)
draw(h,column_title="LINEs")

h <- plot_heat(filter(tes,dfam_name %in% rownames(sce2) & classification  %in% c("LINE","LTR"))$dfam_name,order_by_cl3 = T)
draw(h,column_title="LINEs & LTRs")

h <- plot_heat(filter(tes,dfam_name %in% rownames(sce2) & classification  %in% c("LINE","LTR"))$dfam_name)
draw(h,column_title="LINEs & LTRs")



dev.off()

writexl::write_xlsx(as.data.frame(assay(sce2,"scaled")),snakemake@output$xlsx)
