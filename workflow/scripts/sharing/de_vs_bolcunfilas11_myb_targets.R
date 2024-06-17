Sys.setenv(R_PROFILE=".Rprofile")
source(Sys.getenv("R_PROFILE"))

library(tidyverse)
library(scater)
library(scran)
library(scuttle)
library(fgsea)
library(tesseract)
library(pdftools)

pngfile <- pdftools::pdf_convert("data/bolcun-filas2011_dev067645tables3.pdf",dpi=600)

text1 <- tesseract::ocr(pngfile[1])
text2 <- tesseract::ocr(pngfile[2])
text3 <- tesseract::ocr(pngfile[3])


get_genes <- \(x) str_split(x,"\n(?=chr)") |> paste(collapse = "\n") |>
  str_extract_all("(?<=(\\+|-|NA) )\\d+\\s\\S+") |>
  unlist() |>
  str_extract_all("(?<=\\d\\s).+") |>
  unlist() |>
  str_remove_all("\\/|gene\\\\")
  
  
myb_targets <- list(text1,text2,text3) |> map(get_genes)  |> Reduce(`c`,x=_)


sce <- read_rds(ifelse(exists("snakemake"),
                       snakemake@input$sce,
                       "results/germ_cells/adult.sce.germ_cell.both_genotypes.subclustered.reintegrated.rds"))

lkup <- as_tibble(rowData(sce)[c("gene_name","gene_id")]) |>
  dplyr::rename(feature="gene_id")


dat <- read_tsv(ifelse(exists("snakemake"),
                       snakemake@input$tsv,
                       "results/differential_expression/adult.tbl.germ_cell.reprocessed.de.tsv.gz")) |>
  left_join(lkup)



dat |>
  dplyr::select(celltype,feature,gene_name,logFC) |>
  mutate(myb.target = gene_name %in% myb_targets) |>
  ggplot(aes(myb.target,logFC)) +
  geom_boxplot() +
  facet_wrap(~celltype) +
  ggpubr::stat_compare_means()

ranks <- split(dat,dat$celltype) |> map(~{dplyr::select(.x,gene_name,logFC)}) |> map(deframe)

gs <- list(bolcun_filas=myb_targets)
gsea_res <- ranks |>
  map(~fgsea(gs,stats=.x,nproc=1,eps=0))

gsea_res |>
  map_df(as_tibble,.id="label") |>
  mutate(padj = p.adjust(pval,method="BH"))  |>
  dplyr::select(label,pathway,pval,padj,NES) |>
  gt::gt()
