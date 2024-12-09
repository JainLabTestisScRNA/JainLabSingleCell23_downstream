Sys.setenv(R_PROFILE=".Rprofile")
source(Sys.getenv("R_PROFILE"))

library(tidyverse)
library(scater)
library(scran)
library(scuttle)
library(ggpubr)

lines <- read_tsv(ifelse(exists("snakemake"),snakemake@input$classification,"data/dfam_classification.parsed.txt")) |>
  filter(classification == "LINE")

sce <- read_rds(ifelse(exists("snakemake"),snakemake@input$sce,"results/germ_cells/adult.sce.germ_cell.both_genotypes.subclustered.reintegrated.rds"))

sce <- sce[rownames(sce) %in% lines$dfam_name,]
  
sce2 <- aggregateAcrossFeatures(sce,ids=rep("LINE",nrow(sce)))

sce2 <- logNormCounts(sce2)

g <- makePerCellDF(sce2,features="LINE") |>
  filter(label!="8/Elongating") |>
  ggplot(aes(genotype,LINE,fill=genotype)) +
  #geom_boxplot() +
  geom_violin(scale = "width") +
  stat_summary(geom="point",color="red",position = position_dodge(width=0.8)) +
  ylab("LINE expression\n(summed + normalized)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle=45,hjust=1)) +
  guides(fill="none") +
  scale_fill_grey() +
  ggpubr::stat_compare_means(ref.group = "WT",label = "p.format",size=2.5) +
  facet_wrap(~label,nrow = 1)

ggsave(snakemake@output$pdf, g, height = 3.5,width = 9)

write_tsv(g$data,snakemake@output$tsv)



#makePerCellDF(sce,features="L1MdA_I_5end") |>
#  ggplot(aes(L1MdA_I_5end,fill=genotype)) +
#  geom_histogram(position = "stack") +
#  facet_grid(genotype~label,scales="free") +
#  scale_y_log10()
