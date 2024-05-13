Sys.setenv(R_PROFILE=".Rprofile")
source(Sys.getenv("R_PROFILE"))

library(tidyverse)
library(scuttle)
library(scran)
library(SingleCellExperiment)

theme_set(theme_classic())

te_classes <- ifelse(exists("snakemake"),
                     snakemake@input$cl,
                     "data/dfam_classification.parsed.txt")

te_classes <- te_classes |> read_tsv() |> dplyr::select(feature="dfam_name",classification)


sce_fl <- ifelse(exists("snakemake"),
             snakemake@input$sce,
             "results/differential_expression/adult.sce.germ_cell.reprocessed.pseudobulk.rds")

sce <- read_rds(sce_fl)

lkup <- rowData(sce)[,c("gene_name"),drop=F] |>
  as_tibble(rownames="feature")


fl <- ifelse(exists("snakemake"),
             snakemake@input$tsv,
             "results/differential_expression/adult.tbl.germ_cell.reprocessed.de.tsv.gz")

de <- read_tsv(fl)

de <- de |> left_join(lkup, by="feature") |>
  dplyr::relocate(gene_name,.after=feature)

de<- de |> 
  mutate(feat.type = if_else(str_detect(feature,"ENSMUSG"),"gene","TE"))


de <- de |> left_join(te_classes,by="feature") |> mutate(classification =replace_na(classification,"gene"))

# handle ordering whther using broad anno or numbered germ clusters
if ("Spermatogonia" %in% de$celltype) {
  celltype_ord <- c("Spermatogonia","Spermatocyte","RoundSpermatid","Elongating") 
} else {
  celltype_ord <- unique(as.character(de$celltype))
  celltype_ord <- celltype_ord[str_extract(celltype_ord,"\\d+") |> as.integer() |> order()]
}

de <- de |>
  mutate(celltype = fct_relevel(as.character(celltype),celltype_ord))

logcpmlims <- c(min(de$logCPM),max(de$logCPM))
#logfclims <- c(-max(abs(de$logFC)),max(abs(de$logFC)))
logfclims <- c(-3,3)
logfdrlims <- c(0,7)


plot_ma <- function(x) {
  arrange(x,-str_detect(feature,"ENSMUSG")) |>
    mutate(oob=abs(logCPM) > logcpmlims[2] | abs(logFC) > logfclims[2]) |>
    ggplot(aes(logCPM,logFC)) +
    geom_point(data=\(x)filter(x,feat.type!="TE"&FDR >0.05 &!oob),color="gray",size=rel(0.25)) +
    geom_point(data=\(x)filter(x,feat.type!="TE"&FDR <=0.05&!oob),color="gray48",size=rel(0.25)) +
    geom_point(data=\(x)filter(x,classification=="LINE"&FDR >0.05&!oob),color="lightcoral",size=rel(0.25)) +    
    geom_point(data=\(x)filter(x,classification=="LINE"&FDR <=0.05&!oob),color="red",size=rel(0.5)) +
    geom_point(data=\(x)filter(x,classification=="LTR"&FDR >0.05&!oob),color="royalblue1",size=rel(0.25)) +    
    geom_point(data=\(x)filter(x,classification=="LTR"&FDR <=0.05&!oob),color="blue",size=rel(0.5)) +
    geom_point(data=\(x)filter(x,classification=="DNA"&FDR >0.05&!oob),color="palegreen1",size=rel(0.25)) +    
    geom_point(data=\(x)filter(x,classification=="DNA"&FDR <=0.05&!oob),color="darkgreen",size=rel(0.5)) +
    geom_point(data=\(x)filter(x,classification=="SINE"&FDR >0.05&!oob),color="goldenrod1",size=rel(0.25)) +    
    geom_point(data=\(x)filter(x,classification=="SINE"&FDR <=0.05&!oob),color="darkorange3",size=rel(0.5)) +
    
    geom_point(data=\(x)filter(x,feat.type!="TE"&FDR >0.05 &oob),color="black",size=rel(1.25),shape=15) +
    geom_point(data=\(x)filter(x,feat.type!="TE"&FDR <=0.05&oob),color="black",size=rel(1.25),shape=15) +
    geom_point(data=\(x)filter(x,classification=="LINE"&FDR >0.05&oob),color="lightcoral",size=rel(1.25),shape=15) +    
    geom_point(data=\(x)filter(x,classification=="LINE"&FDR <=0.05&oob),color="red",size=rel(1.25),shape=15) +
    geom_point(data=\(x)filter(x,classification=="LTR"&FDR >0.05&oob),color="royalblue1",size=rel(1.25),shape=15) +    
    geom_point(data=\(x)filter(x,classification=="LTR"&FDR <=0.05&oob),color="blue",size=rel(1.25),shape=15) +
    geom_point(data=\(x)filter(x,classification=="DNA"&FDR >0.05&oob),color="palegreen1",size=rel(1.25),shape=15) +    
    geom_point(data=\(x)filter(x,classification=="DNA"&FDR <=0.05&oob),color="darkgreen",size=rel(1.25),shape=15) +
    geom_point(data=\(x)filter(x,classification=="SINE"&FDR >0.05&oob),color="goldenrod1",size=rel(1.25),shape=15) +    
    geom_point(data=\(x)filter(x,classification=="SINE"&FDR <=0.05&oob),color="darkorange3",size=rel(1.25),shape=15) +
    
    ggrepel::geom_text_repel(data = \(x) filter(x,str_detect(feature,"L1MdA_I_")),
                             aes(label=feature),seed = 2,
                             force = 1000,
                             nudge_y = 0.5,nudge_x = 0.5,
                             box.padding = 2,
                             min.segment.length = unit(0, 'lines'),
                             color="red",fontface="bold") +
    facet_wrap(~celltype,ncol=4) +
    geom_hline(yintercept=0,linetype="dotted") +
    theme(aspect.ratio = 0.5) +
    ylab("logFC (MUT/WT)") +
    #coord_cartesian(xlim=logcpmlims,ylim=logfclims) +
    scale_y_continuous(oob=scales::squish,limits = logfclims) +
    scale_x_continuous(oob=scales::squish,limits = logcpmlims) 
}

plot_ma <- function(x) {
  arrange(x,-str_detect(feature,"ENSMUSG")) |>
    mutate(oob=abs(logCPM) > logcpmlims[2] | abs(logFC) > logfclims[2]) |>
    ggplot(aes(logCPM,logFC)) +
    geom_point(data=\(x)filter(x,feat.type!="TE"&FDR >0.05 &!oob),color="gray",size=rel(0.25)) +
    geom_point(data=\(x)filter(x,feat.type!="TE"&FDR <=0.05&!oob),color="gray48",size=rel(0.25)) +
    geom_point(data=\(x)filter(x,classification=="LINE"&FDR >0.05&!oob),color="lightcoral",size=rel(0.25)) +    
    geom_point(data=\(x)filter(x,classification=="LINE"&FDR <=0.05&!oob),aes(color="LINE"),size=rel(0.5)) +
    geom_point(data=\(x)filter(x,classification=="LTR"&FDR >0.05&!oob),color="royalblue1",size=rel(0.25)) +    
    geom_point(data=\(x)filter(x,classification=="LTR"&FDR <=0.05&!oob),aes(color="LTR"),size=rel(0.5)) +
    geom_point(data=\(x)filter(x,classification=="DNA"&FDR >0.05&!oob),color="palegreen1",size=rel(0.25)) +    
    geom_point(data=\(x)filter(x,classification=="DNA"&FDR <=0.05&!oob),aes(color="DNA"),size=rel(0.5)) +
    geom_point(data=\(x)filter(x,classification=="SINE"&FDR >0.05&!oob),color="goldenrod1",size=rel(0.25)) +    
    geom_point(data=\(x)filter(x,classification=="SINE"&FDR <=0.05&!oob),aes(color="SINE"),size=rel(0.5)) +
    
    geom_point(data=\(x)filter(x,feat.type!="TE"&FDR >0.05 &oob),color="black",size=rel(1.25),shape=15) +
    geom_point(data=\(x)filter(x,feat.type!="TE"&FDR <=0.05&oob),color="black",size=rel(1.25),shape=15) +
    geom_point(data=\(x)filter(x,classification=="LINE"&FDR >0.05&oob),color="lightcoral",size=rel(1.25),shape=15) +    
    geom_point(data=\(x)filter(x,classification=="LINE"&FDR <=0.05&oob),aes(color="LINE"),size=rel(1.25),shape=15) +
    geom_point(data=\(x)filter(x,classification=="LTR"&FDR >0.05&oob),color="royalblue1",size=rel(1.25),shape=15) +    
    geom_point(data=\(x)filter(x,classification=="LTR"&FDR <=0.05&oob),aes(color="LTR"),size=rel(1.25),shape=15) +
    geom_point(data=\(x)filter(x,classification=="DNA"&FDR >0.05&oob),color="palegreen1",size=rel(1.25),shape=15) +    
    geom_point(data=\(x)filter(x,classification=="DNA"&FDR <=0.05&oob),aes(color="DNA"),size=rel(1.25),shape=15) +
    geom_point(data=\(x)filter(x,classification=="SINE"&FDR >0.05&oob),color="goldenrod1",size=rel(1.25),shape=15) +    
    geom_point(data=\(x)filter(x,classification=="SINE"&FDR <=0.05&oob),aes(color="SINE"),size=rel(1.25),shape=15) +
    
    ggrepel::geom_text_repel(data = \(x) filter(x,str_detect(feature,"L1MdA_I_")),
                             aes(label=feature),seed = 2,
                             force = 1000,
                             nudge_y = 0.5,nudge_x = 0.5,
                             box.padding = 2,
                             min.segment.length = unit(0, 'lines'),
                             color="red",fontface="bold") +
    facet_wrap(~celltype,ncol=4) +
    geom_hline(yintercept=0,linetype="dotted") +
    ylab("logFC (MUT/WT)") +
    #coord_cartesian(xlim=logcpmlims,ylim=logfclims) +
    scale_y_continuous(oob=scales::squish,limits = logfclims) +
    scale_x_continuous(oob=scales::squish,limits = logcpmlims) +
    scale_color_manual(values =c(DNA="darkgreen",LTR="blue",LINE="red",SINE="orange"),name="") +
    theme(aspect.ratio = 0.5)
}



pdf(snakemake@output$ma,height = 3,width = 5)

split(de,de$celltype) |>
  map(plot_ma)
  
dev.off()

plot_volc <-  function(x) {
    arrange(x,-str_detect(feature,"ENSMUSG")) |>
    mutate(oob=abs(-log10(FDR)) > logfdrlims[2] | abs(logFC) > logfclims[2]) |>
    ggplot(aes(logFC,-log10(FDR))) +
    geom_point(data=\(x)filter(x,feat.type!="TE"&FDR >0.05 &!oob),color="gray",size=rel(0.25)) +
    geom_point(data=\(x)filter(x,feat.type!="TE"&FDR <=0.05&!oob),color="gray48",size=rel(0.25)) +
    geom_point(data=\(x)filter(x,feat.type=="TE"&FDR >0.05&!oob),color="firebrick1",size=rel(0.25)) +    
    geom_point(data=\(x)filter(x,feat.type=="TE"&FDR <=0.05&!oob),color="red",size=rel(0.5)) +
    
    geom_point(data=\(x)filter(x,feat.type!="TE"&FDR >0.05 &oob),color="gray",size=rel(1.25),shape=15) +
    geom_point(data=\(x)filter(x,feat.type!="TE"&FDR <=0.05&oob),color="gray48",size=rel(1.25),shape=15) +
    geom_point(data=\(x)filter(x,feat.type=="TE"&FDR >0.05&oob),color="firebrick1",size=rel(1.25),shape=15) +    
    geom_point(data=\(x)filter(x,feat.type=="TE"&FDR <=0.05&oob),color="red",size=rel(1.25),shape=15) +
    
    ggrepel::geom_text_repel(data = \(x) filter(x,str_detect(feature,"L1MdA_I_")),seed = 2,
                             aes(label=feature),
                             force = 100,
                             nudge_y = 0.5,nudge_x = 0.5,
                             box.padding = 1,
                             min.segment.length = unit(0, 'lines'),
                             color="red",fontface="bold") +
    facet_wrap(~celltype,ncol=1) +
    geom_hline(yintercept=-log10(0.05),linetype="dotted") +
    xlab("logFC (MUT/WT)") +
    scale_y_continuous(oob=scales::squish,limits = logfdrlims) +
    scale_x_continuous(oob=scales::squish,limits = logfclims) 
}

pdf(snakemake@output$volcano,width = 4,height = 4)

split(de,de$celltype) |>
  map(plot_volc)

dev.off()

write_tsv(de,snakemake@output$tsv)
