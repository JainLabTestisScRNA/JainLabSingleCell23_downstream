library(tidyverse)
library(scater)
library(scran)
library(scuttle)

sce <- read_rds(ifelse(exists("snakemake"),
                 snakemake@input$sce,
                 "results/germ_cells/adult.sce.germ_cell.both_genotypes.subclustered.reintegrated.rds"))

sce <- sce[,sce$genotype == "WT"]

sce2 <- aggregateAcrossCells(sce,sce$label,coldata.merge=list(pseudotime=mean),statistics="sum")
sizeFactors(sce2) <- librarySizeFactors(sce2)
sce2 <- logNormCounts(sce2)

tes <- read_tsv(ifelse(exists("snakemake"),snakemake@input$classifications,"data/dfam_classification.parsed.txt"))# |>
  #filter(classification == "LINE")

te_lkup <- tes |> dplyr::select(dfam_name,classification) |> deframe()
  
sce2 <- sce2[rownames(sce2)%in%tes$dfam_name,]

sce2 <- sce2[rowSums(counts(sce2)) > 100,]

rowData(sce2)$type <- map_chr(rownames(sce2),~{te_lkup[.x]})

pdf(snakemake@output$pdf)
g <- plotHeatmap(sce2,
            features=rownames(sce2),
            center = T,
            scale = T,border_color=NA,
            order_columns_by = c("label","pseudotime"),
            zlim = c(-3,3),
            color_rows_by = "type",
            show_rownames=F)

g
dev.off()

logcounts(sce2) |> t() |> scale() |> t() |>
  as_tibble(rownames="feature") |>
  write_tsv(snakemake@output$tsv)
