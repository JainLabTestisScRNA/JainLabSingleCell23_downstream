library(tidyverse)
library(scater)
library(scran)
library(scuttle)

theme_set(theme_classic())
 
de_by <-  ifelse(exists("snakemake"),
                 snakemake@params$de_by,
                 "celltype")
# regex for column used for differential expression (currently 'label' for GC clusters,"celltype" for overall)
de_by <- sprintf("^%s$",de_by)

sce_fl <- ifelse(exists("snakemake"),
             snakemake@input$sce,
             "results/integration/adult.sce.integrated.clustered.celltypes.rds")

sce <- read_rds(sce_fl)

de_fl <- ifelse(exists("snakemake"),
                snakemake@input$de,
                "results/differential_expression/adult.tbl.pseudobulk.de.tsv.gz")

de <- read_tsv(de_fl)

classification_fl <- 
  ifelse(exists("snakemake"),
         snakemake@input$classification,
         "data/dfam_classification.parsed.txt")

classification <- read_tsv(classification_fl) |>
  dplyr::select(dfam_name,classification)

df <- filter(de,!str_detect(feature,"ENSMUSG") & FDR < 0.05 & logFC > 0) |>
  left_join(classification,by=c(feature="dfam_name")) |>
  dplyr::select(celltype,feature,classification)

df2 <- makePerCellDF(sce,
                     features = df$feature,assay.type = "logcounts") |>
  as_tibble(rownames="cell") |>
  dplyr::select(cell,genotype,batch,celltype=matches(de_by),Barcode,df$feature) |>
  pivot_longer(-c(cell,Barcode,genotype,batch,celltype),names_to = "feature",values_to = "logcounts")

df3 <- inner_join(df2,df,by = join_by(celltype,feature))

df4 <- df3 |> 
  arrange(celltype) |>
  mutate(split_group = paste(celltype,":",feature))

# per DJ's recommendation
df4$simple_batch <- df4$batch |> str_remove("_\\d$")

retransformed_untransformed_mean <- \(x) log2(mean((2^x)-1)+1)

plot_jitters <- function(x) {
  title <- unique(x$split_group)
  ggplot(x,aes(genotype,logcounts,color=genotype)) +
    geom_point(position=position_jitterdodge(seed=2),size=0.2) +
    geom_violin(scale="width",fill=NA) +
    facet_wrap(~feature + simple_batch,scales = "free",nrow=1) +
    stat_summary(geom="point",fun = mean,
                 color="red",aes(group=genotype),position=position_dodge(width = 1)) +
    scale_color_grey(start = 0,end=0.4) +
    guides(color="none")  +
    labs(title=title) +
    theme(aspect.ratio = 0.5)
}

pdf(snakemake@output$pdf)

split(df4,df4$split_group) |>
  map(plot_jitters)

dev.off()

write_tsv(df4,snakemake@output$tsv)

