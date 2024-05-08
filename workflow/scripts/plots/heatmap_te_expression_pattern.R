Sys.setenv(R_PROFILE=".Rprofile")
source(Sys.getenv("R_PROFILE"))

library(tidyverse)
library(scater)
library(scran)
library(scuttle)
library(ComplexHeatmap)
library(RColorBrewer)

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
sce2 <- sce2[rowSums(counts(sce2)) > 10,]

# include info about te types
rowData(sce2)$type <- map_chr(rownames(sce2),~{te_lkup[.x]})
rd <- rowData(sce2)[,"type",drop=F] |> as.data.frame()

assay(sce2,"scaled") <- logcounts(sce2) |> t() |> scale() |> t()


mats <- list(WT=assay(sce2[,sce2$genotype == "WT"],"scaled"),
             MUT=assay(sce2[,sce2$genotype == "MUT"],"scaled"))

plot_heat <-  function(gt,feats) {
  # row ord/clustering based on both gt
  rdend <- hclust(dist(assay(sce2,"scaled")[feats,]))
  m <- mats[[gt]][feats,]
  m <- m[rdend$order,]
  Heatmap(m,cluster_rows = F,cluster_columns=F,show_row_names = F,show_column_names = T,col = rev(coul),name = "Expression\n(row z-scores)",
          row_dend_width = unit(20, "mm"), split=factor(rd[rownames(m),"type"]),border = "black",
          column_title=gt,column_labels = colnames(m))
}

plot_heat("WT",filter(tes,dfam_name %in% rownames(sce2) & classification  %in% c("LINE","LTR"))$dfam_name) +
  plot_heat("MUT",filter(tes,dfam_name %in% rownames(sce2) & classification  %in% c("LINE","LTR"))$dfam_name)

pdf(snakemake@output$pdf)

plot_heat("WT") + plot_heat("MUT")

plot_heat("WT",filter(tes,dfam_name %in% rownames(sce2) & classification  %in% c("LINE","LTR"))$dfam_name) +
  plot_heat("MUT",filter(tes,dfam_name %in% rownames(sce2) & classification  %in% c("LINE","LTR"))$dfam_name)




dev.off()


writexl::write_xlsx(map(mats,as.data.frame),snakemake@output$xlsx)
