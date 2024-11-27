Sys.setenv(R_PROFILE=".Rprofile")
source(Sys.getenv("R_PROFILE"))

library(tidyverse)
library(scater)
library(scran)
library(scuttle)

theme_set(theme_classic())

fl <- ifelse(exists("snakemake"),snakemake@input$sce,"results/integration/adult.sce.integrated.clustered.celltypes.rds")
sce <- read_rds(fl)

marker.info <- scoreMarkers(sce, sce$celltype)

chosen <- marker.info$InnateLymphoid

ordered <- chosen[order(chosen$mean.AUC, decreasing=TRUE),]
ordered$gene_name <- rowData(sce)["gene_name"][rownames(ordered),]

head(ordered[,c("gene_name","mean.AUC")],20) |> as_tibble() |> print(n=Inf)


g <- plotExpression(sce,x="celltype",other_fields = "celltype",
                    features = c("Id2","Il7r","Rora","Thy1","Ccl5","Cd52"),
                    swap_rownames = "gene_name",scales = "free",one_facet = F) +
  theme(axis.text.x=element_text(angle=45,hjust=1)) +
  facet_wrap(~Feature,ncol=1,scales="free")


g



ggsave(snakemake@output$pdf,g,height = 12)

g$data |> write_tsv(snakemake@output$tsv)
