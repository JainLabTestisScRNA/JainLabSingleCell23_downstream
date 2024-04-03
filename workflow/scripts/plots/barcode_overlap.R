library(scater)
library(scran)
library(scuttle)
library(tidyverse)
library(ggVennDiagram)

#fls <- Sys.glob("results/cellranger-import/adult_*.sce.assigned.rds")
fls <- snakemake@input$fls

names(fls) <- str_extract(fls,"(?<=adult_).+(?=\\.sce)")

bcl <- map(fls,read_rds) |>
  map(`$`,"Barcode")

# ggVennDiagram uses this under the hood
d <- process_region_data(Venn(bcl)) |>
  dplyr::select(-item)

g <- ggVennDiagram(bcl,) + 
  labs(title="barcode overlap among batches (non-empty cells assigned CMO)") + 
  theme(plot.title=element_text(hjust=0.5,face="bold")) +
  scale_fill_gradient(low="white",high="white") +
  scale_color_manual(values = rep("black",length(bcl))) +
  guides(fill="none")


write_tsv(d,snakemake@output$tsv)
ggsave(snakemake@output$pdf,g)
