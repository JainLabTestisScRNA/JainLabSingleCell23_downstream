library(tidyverse)

theme_set(theme_classic())

fl <- ifelse(exists("snakemake"),
             snakemake@input$tsv,
             "results/differential_expression/adult.tbl.germ_cell.reprocessed.de.tsv.gz")

de <- read_tsv(fl)

de<- de |> 
  mutate(feat.type = if_else(str_detect(feature,"ENSMUSG"),"gene","TE"))


# handle ordering whther using broad anno or numbered germ clusters
if ("Spermatogonia" %in% de$celltype) {
  celltype_ord <- c("Spermatogonia","Spermatocyte","RoundSpermatid","Elongating") 
} else {
  celltype_ord <- unique(as.character(de$celltype))
  celltype_ord <- celltype_ord[str_extract(celltype_ord,"\\d+") |> as.integer() |> order()]
}

de <- de |>
  mutate(celltype = fct_relevel(as.character(celltype),celltype_ord))

plot_ma <- function(x) {
  arrange(x,-str_detect(feature,"ENSMUSG")) |>
    ggplot(aes(logCPM,logFC)) +
    geom_point(data=\(x)filter(x,FDR>0.05),color="lightgray",size=rel(0.5)) +
    geom_point(data=\(x)filter(x,FDR<=0.05),color="darkgray",size=rel(0.5)) +
    geom_point(data=\(x)filter(x,feat.type=="TE"&FDR <0.05),size=rel(2),shape=21,color="red",fill=NA) +
    ggrepel::geom_text_repel(data = \(x) filter(x,str_detect(feature,"L1MdA_I_")),aes(label=feature),min.segment.length = 0.1) +
    facet_wrap(~celltype,ncol=4) +
    geom_hline(yintercept=0,linetype="dotted") +
    theme(aspect.ratio = 0.5) +
    ylab("logFC (MUT/WT)")
}



pdf(snakemake@output$ma)

split(de,de$celltype) |>
  map(plot_ma)
  
dev.off()

plot_volc <-  function(x) {
    arrange(x,-str_detect(feature,"ENSMUSG")) |>
    ggplot(aes(logFC,-log10(FDR),color=feat.type)) +
    geom_point(data=\(x)filter(x,FDR>0.05),color="lightgray",size=rel(0.5)) +
    geom_point(data=\(x)filter(x,FDR<=0.05),color="darkgray",size=rel(0.5)) +
    geom_point(data=\(x)filter(x,feat.type=="TE"&FDR <0.05),size=rel(2),shape=21,color="red",fill=NA) +
    facet_wrap(~celltype,ncol=1) +
    geom_hline(yintercept=-log10(0.05),linetype="dotted") +
    theme(aspect.ratio = 0.5) +
    scale_color_manual(values = c(gene="darkgray","TE"="red")) +
    xlab("logFC (MUT/WT)")
}

pdf(snakemake@output$volcano)

split(de,de$celltype) |>
  map(plot_volc)

dev.off()

write_tsv(de,snakemake@output$tsv)