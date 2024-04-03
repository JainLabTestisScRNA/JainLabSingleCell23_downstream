library(tidyverse)
library(scater)
library(scran)
library(scuttle)
library(BiocParallel)
library(tictoc)

threads <- 6
threads <-snakemake@threads
bpp <- MulticoreParam(threads)

#sce_fl <- "results/germ_cells/adult.sce.integrated.clustered.celltypes.germ_cell.reprocessed.rds"
sce_fl <-  snakemake@input$sce
sce <- read_rds(sce_fl)

# see https://bioconductor.org/books/3.18/OSCA.advanced/cell-cycle-assignment.html
mm.pairs <- readRDS(system.file("exdata", "mouse_cycle_markers.rds", 
                                package="scran"))

#sce2 <- sce[,sample(1:ncol(sce),size=1000)]
tic()
cc <- cyclone(sce,mm.pairs,BPPARAM=bpp,assay.type="logcounts",iter=1000)
toc()

sce$cell_cycle <- cc$phases

write_rds(sce,snakemake@output$rds)
write_rds(cc,snakemake@output$cc)
#plotReducedDim(sce,dimred="corrected",ncomponents = c(1,3),colour_by = "cell_cycle") + facet_wrap(~colour_by)
