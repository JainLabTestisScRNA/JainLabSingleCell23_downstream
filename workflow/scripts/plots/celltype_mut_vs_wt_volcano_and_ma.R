library(tidyverse)

theme_set(theme_classic())

germ_cell_types <- c("Spermatogonia","Spermatocyte","RoundSpermatid","Elongating")

fl <- ifelse(exists("snakemake"),
             snakemake@input$tsv,
             "results/differential_expression/adult.tbl.pseudobulk.de.tsv.gz")

de <- read_tsv(fl)

de.germ <- de |> 
  filter(celltype %in% germ_cell_types) |>
  mutate(feat.type = if_else(str_detect(feature,"ENSMUSG"),"gene","TE"))

pdf(snakemake@output$ma)
de.germ |>
  arrange(-str_detect(feature,"ENSMUSG")) |>
  ggplot(aes(logCPM,logFC,color=feat.type)) +
  geom_point(size=rel(0.5)) +
  geom_point(data=\(x)filter(x,feat.type=="TE"&FDR <0.05),size=rel(2),shape=21,color="black",fill=NA) +
  facet_wrap(~celltype,ncol=1) +
  geom_hline(yintercept=0,linetype="dotted") +
  theme(aspect.ratio = 0.5) +
  scale_color_manual(values = c(gene="gray","TE"="red")) +
  ylab("logFC (MUT/WT)")
dev.off()

pdf(snakemake@output$volcano)
de.germ |>
  arrange(-str_detect(feature,"ENSMUSG")) |>
  ggplot(aes(logFC,-log10(FDR),color=feat.type)) +
  geom_point(size=rel(0.5)) +
  geom_point(data=\(x)filter(x,feat.type=="TE"&FDR <0.05),size=rel(2),shape=21,color="black",fill=NA) +
  facet_wrap(~celltype,ncol=1) +
  geom_hline(yintercept=-log10(0.05),linetype="dotted") +
  theme(aspect.ratio = 0.5) +
  scale_color_manual(values = c(gene="gray","TE"="red")) +
  xlab("logFC (MUT/WT)")
dev.off()

write_tsv(de.germ,snakemake@output$tsv)