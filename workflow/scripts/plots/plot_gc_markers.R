Sys.setenv(R_PROFILE=".Rprofile")
source(Sys.getenv("R_PROFILE"))

library(tidyverse)
library(scater)
library(scran)
library(scuttle)

sce <- read_rds(ifelse(exists("snakemake"),snakemake@input$sce,"results/germ_cells/adult.sce.germ_cell.wt.subclustered.reintegrated.rds"))
markers <- read_tsv(ifelse(exists("snakemake"),snakemake@input$tsv,"results/germ_cells/adult.sce.germ_cell.wt.subclustered.reintegrated.markers.tsv.gz"))

## from osca
## Cohen’s  d is a standardized log-fold change where the difference in the mean
## log-expression between groups is scaled by the average standard deviation across
## groups. In other words, it is the number of standard deviations that separate 
## the means of the two groups. The interpretation is similar to the log-fold change; 
## positive values indicate that the gene is upregulated in our cluster of interest, 
## negative values indicate downregulation and values close to zero indicate that 
## there is little difference. Cohen’s d is roughly analogous to the t-statistic 
## in various two-sample t-tests.
## ...
## The minimum rank, a.k.a., “min-rank” (rank.*) is the smallest rank of each gene across all pairwise comparisons. 

chosen <- markers |>
  filter(str_detect(feature,"ENSMUSG")) |>
  group_by(label) |>
  slice_min(rank.logFC.cohen,n=20) |>
  slice_max(min.logFC.cohen,n=20) |>
  ungroup()

pdf(snakemake@output$pdf,width=8.5,height=11)

plotGroupedHeatmap(sce, features=unique(chosen$gene_name), group="label", block = "batch",scale = T,
                   center=TRUE, zlim=c(-3, 3),swap_rownames = "gene_name",cluster_cols=F,border_color=NA,show_rownames=T)

dev.off()

write_tsv(chosen,snakemake@output$tsv)