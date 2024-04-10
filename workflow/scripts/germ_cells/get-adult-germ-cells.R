library(tidyverse)
library(scater)
library(scran)
library(scuttle)

fl <- "results/integration/adult.sce.integrated.clustered.celltypes.rds"
fl <- snakemake@input$sce

sce <- read_rds(fl)

# these data will require reidentifying highly variable genes and pca
assay(sce,"reconstructed") <- NULL
reducedDims(sce) <- NULL
sce$label2 <- NULL
sce$celltype.garnett <- NULL
rowSubset(sce) <- NULL

# see methods of green et al 'Focused analysis for germ cells'
sce <- sce[,sce$nexprs >= 1000]

# subset data - all germ cells
germ_types <- c("Spermatogonia","Spermatocyte","RoundSpermatid","Elongating")
sce2 <- sce[,sce$celltype %in% germ_types]
write_rds(sce2,snakemake@output$rds)

# subset data all mut germs cells
sce2.mut <- sce2[,sce2$genotype == "MUT"]
write_rds(sce2.mut,snakemake@output$mut_rds)

# subset data all wt germs cells
sce2.wt <- sce2[,sce2$genotype == "WT"]
write_rds(sce2.wt,snakemake@output$wt_rds)


sce2.wt.spg <- sce2[,sce2$genotype == "WT" & sce2$celltype == "Spermatogonia"]
write_rds(sce2.wt.spg,snakemake@output$wt_spg_rds)

sce2.wt.spc <- sce2[,sce2$genotype == "WT" & sce2$celltype == "Spermatocyte"]
write_rds(sce2.wt.spc,snakemake@output$wt_spc_rds)

sce2.wt.rspt <- sce2[,sce2$genotype == "WT" & sce2$celltype == "RoundSpermatid"]
write_rds(sce2.wt.rspt,snakemake@output$wt_rspt_rds)

sce2.wt.espt <- sce2[,sce2$genotype == "WT" & sce2$celltype == "Elongating"]
write_rds(sce2.wt.espt,snakemake@output$wt_espt_rds)






