Sys.setenv(R_PROFILE=".Rprofile")
source(Sys.getenv("R_PROFILE"))

library(tidyverse)
library(scater)
library(scran)
library(scuttle)
library(fgsea)
library(GSEABase)
library(writexl)

gs <- read_rds(ifelse(exists("snakemake"),snakemake@input$gsea,"results/enrichment/gsea_prerank.rds"))
ranks <- read_rds(ifelse(exists("snakemake"),snakemake@input$ranks,"results/enrichment/gsea_prerank.ranks.rds"))
gmt <- gmtPathways(ifelse(exists("snakemake"),snakemake@input$msigdb,"data/m5.all.v2023.2.Mm.symbols.gmt"))



plot_table <-  function(comparison,n=10) {
  rnk <- ranks[[comparison]]
  pwPos <- head(gs[[comparison]][sign(NES)==1, ][order(pval),pathway],n)
  pwNeg <- head(gs[[comparison]][sign(NES)==-1, ][order(pval),pathway],n)
  pw<- gmt[c(pwPos,rev(pwNeg))]
  res <- gs[[comparison]] |> 
    as_tibble() |> 
    filter(pathway %in% names(pw)) |> 
    mutate(leadingEdge=map_chr(leadingEdge,paste,collapse=",")) |>
    arrange(-NES)
  g <- fgsea::plotGseaTable(stats=rnk, fgseaRes = gs[[comparison]],pathways = pw) + labs(tag=comparison)
  return(list(data=res,g=g))
}

res <- names(ranks) |>
  set_names() |>
  map(plot_table)

dat <- res |> map(pluck,"data")
names(dat) <- names(dat) |> str_replace("\\/",".")
write_xlsx(dat,snakemake@output$xlsx)

pdf(snakemake@output$pdf)

res |> map(pluck,"g")

dev.off()