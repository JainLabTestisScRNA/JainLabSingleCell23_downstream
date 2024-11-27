Sys.setenv(R_PROFILE=".Rprofile")
source(Sys.getenv("R_PROFILE"))

library(tidyverse)
library(scater)
library(scran)
library(scuttle)

theme_set(theme_classic())

fl <- ifelse(exists("snakemake"),snakemake@input$sce,"results/germ_cells/adult.sce.integrated.clustered.celltypes.germ_cell.reprocessed.rds")

sce <- read_rds(fl)

g <- plotExpression(sce,x = "label",
                    features = c("Sycp3","Stra8","Dmrt1"),
                    swap_rownames = "gene_name",scales = "free",ncol = 1) +
  theme(axis.text.x=element_text(angle=45,hjust=1))

# put x axis into the correct order
labs <- sce$label |> unique() 
level_order <- labs[as.integer(str_extract(labs,"\\d+")) |> order()]
g <- g + scale_x_discrete(limits = level_order)

ggsave(snakemake@output$pdf,g)

g$data |> write_tsv(snakemake@output$tsv)

# ------------
# oi <- list(Sycp3="ENSMUSG00000020059",Stra8="ENSMUSG00000029848",Dmrt1="ENSMUSG00000024837")
# 
# x <- makePerCellDF(sce,features = oi) |>
#   as_tibble(rownames="cell") |>
#   dplyr::select(cell,label,batch,genotype,contains("ENSMUSG"))
# 
# 
# preL <- c("")
# x |>
#   filter(label %in% preL) |>
#   ggplot(aes(ENSMUSG00000020059,ENSMUSG00000029848,color=label)) +
#   geom_point()
# 
# 
# x |>
#   filter(label%in%preL) |>
#   ggplot(aes(ENSMUSG00000020059,ENSMUSG00000024837,color=label)) +
#   geom_point()
# 
# 
# x |>
#   filter(label%in%preL) |>
#   ggplot(aes(ENSMUSG00000029848,ENSMUSG00000024837,color=label)) +
#   geom_point()
