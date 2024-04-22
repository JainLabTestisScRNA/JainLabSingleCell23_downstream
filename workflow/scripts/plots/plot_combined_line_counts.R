library(tidyverse)
library(scater)
library(scran)
library(scuttle)

lines <- read_tsv(ifelse(exists("snakemake"),snakemake@input$classification,"data/dfam_classification.parsed.txt")) |>
  filter(classification == "LINE")

sce <- read_rds(ifelse(exists("snakemake"),snakemake@input$sce,"results/germ_cells/adult.sce.germ_cell.both_genotypes.subclustered.reintegrated.rds"))

sce <- sce[rownames(sce) %in% lines$dfam_name,]
  
sce2 <- aggregateAcrossFeatures(sce,ids=rep("LINE",nrow(sce)))

sce2 <- logNormCounts(sce2)

g <- plotExpression(sce2,features="LINE",x="label",other_fields = "genotype") +
  facet_wrap(~genotype)

ggsave(snakemake@output$pdf, g)

write_tsv(g$data,snakemake@output$tsv)


#makePerCellDF(sce,features="L1MdA_I_5end") |>
#  ggplot(aes(L1MdA_I_5end,fill=genotype)) +
#  geom_histogram(position = "stack") +
#  facet_grid(genotype~label,scales="free") +
#  scale_y_log10()
