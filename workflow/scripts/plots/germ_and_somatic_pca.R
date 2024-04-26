Sys.setenv(R_PROFILE=".Rprofile")
source(Sys.getenv("R_PROFILE"))

library(tidyverse)
library(scuttle)
library(scran)
library(scater)

theme_set(theme_classic())

fl <- ifelse(exists("snakemake"),snakemake@input$sce,"results/integration/adult.sce.integrated.clustered.celltypes.rds")
sce <- read_rds(fl)

odir <- ifelse(exists("snakemake"),snakemake@params$odir,"results/plots/germ_and_somatic_pca/")
dir.create(odir,recursive = T)

germ_types <- c("Spermatogonia","Spermatocyte","RoundSpermatid","Elongating")

x <- makePerCellDF(sce) |> 
  as_tibble(rownames="cell") |> 
  dplyr::select(cell,genotype,celltype,corrected.1,corrected.2,corrected.3) |>
  mutate(celltype_group= if_else(celltype %in% germ_types,"germ","somatic"))

g_pc1_pc2 <- ggplot(arrange(x,celltype_group), aes(corrected.1,corrected.2,color=celltype_group)) +
  geom_point()

g_pc1_pc3 <- ggplot(arrange(x,celltype_group), aes(corrected.1,corrected.3,color=celltype_group)) +
  geom_point()

ggsave(snakemake@output$g_pc1_pc2,g_pc1_pc2)
ggsave(snakemake@output$g_pc1_pc3,g_pc1_pc3)
write_tsv(x,snakemake@output$tsv)

