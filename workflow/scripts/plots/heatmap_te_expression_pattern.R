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

# get wt and mutant mats
mats <- c("WT","MUT") |>
  set_names() |>
  map(~{
    x <- sce2[,sce2$genotype == .x]
    colnames(x) <- make.unique(x$label)
    mat <- logcounts(x)
    mat <- t(scale(t(mat)))
    #mat[order(mat[,"4/Spermatocyte"],decreasing = T),]
    return(mat)
  })


# row ord/clustering based on WT to keep consistent
rdend <- hclust(dist(mats$WT))

plot_heat <-  function(gt) {
  m <- mats[[gt]][rdend$order,]
  Heatmap(m,cluster_rows = F,cluster_columns=F,show_row_names = F,col = rev(coul),name = "Expression\n(row z-scores)",
          row_dend_width = unit(20, "mm"), split=factor(rd[rownames(m),"type"]),border = "black",
          column_title=gt)
}

pdf(snakemake@output$pdf)

plot_heat("WT") + plot_heat("MUT")

dev.off()


writexl::write_xlsx(map(mats,as.data.frame),snakemake@output$xlsx)
