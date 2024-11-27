Sys.setenv(R_PROFILE=".Rprofile")
source(Sys.getenv("R_PROFILE"))

library(tidyverse)
library(scater)
library(scran)
library(scuttle)
library(diptest)

sce <- read_rds(ifelse(exists("snakemake"),snakemake@input$sce,
                       "results/germ_cells/adult.sce.germ_cell.both_genotypes.subclustered.reintegrated.rds"))

#sce <- sce[,sce$label == "4/Spermatocyte"]
sce <- sce[rowSums(counts(sce)) > 100,]

# ------------------------------------------------------------------------------
# proportion 0s
# ------------------------------------------------------------------------------
dat <- counts(sce) |>
  as.matrix() |>
  as_tibble(rownames="feature") |>
  pivot_longer(-feature,names_to = "cell")

dat <- as_tibble(colData(sce)[,c("genotype","label")],rownames="cell") |>
  left_join(dat,y=_)

dat <- dat |>
  group_by(feature,label,genotype) |>
  summarize(prop.zero = sum(value == 0) / n(),.groups = "drop")

gc()

dat <- rowRanges(sce) |>
  as_tibble() |>
  dplyr::select(seqnames,gene_id,gene_name) |>
  left_join(dat,by=c(gene_id = "feature"))


write_rds(dat,snakemake@output$rds)
# 
# # overall look is promising
# dat |>
#   filter(str_detect(label,"Sperm")) |>
#   #filter(!str_detect(seqnames,"chrX|chrY")) |>
#   pivot_wider(names_from = genotype,values_from = "prop.zero") |>
#   ggplot(aes(MUT,WT)) +
#   ggdensity::geom_hdr_points() +
#   geom_abline(slope=1,intercept = 0,linetype="dashed",color="red") +
#   facet_wrap(~label,scales="free")
# 
# # this pretty much mimics the analysis I've shown previously for per-chrom sharing/msci
# # as a sanity check
# dat |>
#   filter(str_detect(label,"Sperm") & str_detect(seqnames,"chrX|chrY|chr9$")) |>
#   pivot_wider(names_from = genotype,values_from = "prop.zero") |>
#   ggplot(aes(MUT,WT)) +
#   ggdensity::geom_hdr_points() +
#   geom_abline(slope=1,intercept = 0,linetype="dashed",color="red") +
#   facet_grid(seqnames~label)
# 
# 
# de <- read_tsv("results/differential_expression/adult.tbl.germ_cell.reprocessed.de.tsv.gz") |> dplyr::rename(label="celltype",gene_id="feature")
# 
# # stuff with a high share score, but positive or low absolute FC (preferably n.s.)
# # could be shared
# # note reasonable concordance with FC for early cells
# dat |>
#   filter(str_detect(label,"Sperm")) |>
#   pivot_wider(names_from = genotype,values_from = "prop.zero") |>
#   inner_join(de) |>
#   mutate(share.score = log2((MUT+0.01)/(WT+0.01))) |>
#   mutate(shared.candidate = (logFC > 0 | FDR > 0.05) & share.score>0.5) |>
#   ggplot(aes(share.score,logFC,color=shared.candidate)) +
#   geom_point(size=rel(0.2)) +
#   facet_wrap(~label,nrow =1)
# 
# dat |>
#   filter(str_detect(label,"Sperm")) |>
#   filter(seqnames %in% c("chrX","chr9")) |>
#   pivot_wider(names_from = genotype,values_from = "prop.zero") |>
#   inner_join(de) |>
#   mutate(share.score = log2((MUT+0.01)/(WT+0.01))) |>
#   mutate(shared.candidate = (logFC > 0 | FDR > 0.05) & share.score>0.5) |>
#   ggplot(aes(seqnames,share.score)) +
#   geom_violin() +
#   stat_summary(geom="point") +
#   facet_wrap(~label,nrow =1) +
#   ggpubr::stat_compare_means(size=3)
# 
# # ------------------------------------------------------------------------------
# # find candidate shared
# # ------------------------------------------------------------------------------
# 
# shared.universe <- dat |>
#   filter(str_detect(label,"Sperm") & str_detect(gene_id,"ENSMUSG")) |>
#   pivot_wider(names_from = genotype,values_from = "prop.zero") |>
#   inner_join(de) |>
#   mutate(share.score = log2((MUT+0.01)/(WT+0.01))) |>
#   mutate(shared.candidate = (logFC > 0 | FDR > 0.05) & share.score>0.5)
# 
# shared <- filter(shared.universe,shared.candidate) |>
#   pull(gene_name) |>
#   unique()
# 
# 
# 
# genoi <- read_csv("data/BhutaniGenoinformativity/abb1723-bhutani-sm-tables-s1-to-s11-zipped/table_s1_mouse_genoinformativity.csv",skip = 1)
# 
# cb_enr <- readxl::read_xlsx("data/meikar2018.xlsx",skip = 2) |>
#   pull(`Gene Symbol`) |>
#   str_extract("\\w+")
# 
# # meikar
# cb_fold_enr <- readxl::read_xlsx("data/meikar2018.xlsx",skip = 2,sheet = 6) |>
#   dplyr::select(gene_name ="name",cb.fold.enr="Fold enrichment")
# 
# # goh et al
# pachy_pirna_targets <- readxl::read_xlsx("~/Downloads/TableS7_628_low-confidence_mRNA_Targets.xlsx")
# pachy_pirna_targets <- readxl::read_xlsx("~/Downloads/TableS1_72_mRNA_Targets.xlsx")
# pachy_pirna_targets <- readxl::read_xlsx("~/Downloads/TableS2_Non-misregulated_mRNAs.xlsx")
# 
# # vourekas et al
# miwi_targets <- readxl::read_xlsx("~/Downloads/41594_2012_BFnsmb2347_MOESM29_ESM.xlsx")
# 
# de |>
#   left_join(as_tibble(rowData(sce)[,c("gene_name","gene_id")])) |>
#   mutate(CB.enriched = gene_name %in% cb_enr,,
#          pachy.targets = gene_name %in% pachy_pirna_targets$Gene,
#          miwi.targets = gene_name %in% miwi_targets$symbol) |>
#   ggplot(aes(CB.enriched,logFC)) +
#   #geom_violin(outlier.shape = NA,scale="width") +
#   geom_boxplot(outlier.shape = NA) +
#   facet_wrap(~label, scales = "free") +
#   ggpubr::stat_compare_means()
# 
# 
# shared.universe |>
#   mutate(CB.enriched = gene_name %in% cb_enr,
#          pachy.targets = gene_name %in% pachy_pirna_targets$Gene,
#          miwi.targets = gene_name %in% miwi_targets$symbol) |>
#   ggplot(aes(CB.enriched,share.score)) +
#   geom_violin(scale = "area") +
#   stat_summary(geom="point",color="red",size=rel(2)) +
#   facet_wrap(~label,scales = "free",nrow=1) +
#   ggpubr::stat_compare_means(size=3)
# 
# dat |>
#   filter(str_detect(label,"Sperm") & str_detect(gene_id,"ENSMUSG")) |>
#   pivot_wider(names_from = genotype,values_from = "prop.zero") |>
#   mutate(share.score = log2((MUT+0.01)/(WT+0.01))) |>
#   inner_join(cb_fold_enr) |>
#   ggplot(aes(cb.fold.enr,share.score)) +
#   geom_point() +
#   facet_wrap(~label) +
#   geom_smooth(method="lm") +
#   ggpubr::stat_cor()
# 
# 
# dplyr::select(genoi,gene_id=Ensembl_gene_ID,Posterior_mean_genoinformativity,GIM_status) |>
#   inner_join(shared.universe) |>
#   ggplot(aes(fct_reorder(GIM_status,logFC),logFC)) +
#   geom_violin() +
#   stat_summary(geom="point") +
#   facet_wrap(~label) +
#   ggpubr::stat_compare_means(label.y = 2)
# 
# dplyr::select(genoi,gene_id=Ensembl_gene_ID,Posterior_mean_genoinformativity,GIM_status) |>
#   inner_join(shared.universe) |>
#   ggplot(aes(fct_reorder(GIM_status,share.score,median),share.score)) +
#   geom_violin() +
#   stat_summary(geom="point") +
#   facet_wrap(~label,nrow=1,scales="free_y") +
#   ggpubr::stat_compare_means(label.y = 2,size=3) +
#   coord_cartesian(ylim=c(-2,2)) +
#   theme(axis.text.x = element_text(angle=45,hjust=1))
# 
# 
# 
# 
# library(fgsea)
# 
# gmt <- gmtPathways(ifelse(exists("snakemake"),snakemake@input$msigdb,"data/m5.all.v2023.2.Mm.symbols.gmt"))
# 
# enr_results <- fora(pathways = gmt,genes = shared,universe = unique(shared.universe$gene_name),minSize = 15,maxSize = 300) |> filter(padj < 0.05) |>
#   as_tibble()
# 
# 
# enr_results |>
#   filter(padj < 0.05) |>
#   slice_min(pval,n=25) |>
#   ggplot(aes(-log10(padj),fct_reorder(pathway,pval))) +
#   geom_col()
# 
# enr_results |> 
#   filter(!str_detect(pathway,"RIBO|TRANSLA|PEPTIDE|BIOSYNTHETIC|METABOLIC|POLYSOME")) |>
#   slice_min(pval,n=20) |>
#   gt::gt()
# 
# 
# 
# 
# 
# # ------------------------------------------------------------------------------
# # dip test
# # ------------------------------------------------------------------------------
# dat <- logcounts(sce) |>
#   as.matrix() |>
#   as_tibble(rownames="feature") |>
#   pivot_longer(-feature,names_to = "cell")
# 
# dat <- as_tibble(colData(sce)[,c("genotype","label")],rownames="cell") |>
#   left_join(dat,y=_)
# 
# #c(rnorm(1000,mean = 10),rnorm(1000,mean = 10)) |> dip.test() |> broom::tidy()
# 
# dat <- dat |>
#   group_by(feature,genotype,label) |>
#   summarise(dip.stat = dip(value,min.is.0 = T),.groups = "drop") |>
#   pivot_wider(names_from = genotype,values_from = "dip.stat") |>
#   mutate(share.score = log2(MUT/WT))
# 
# dat <- rowRanges(sce) |>
#   as_tibble() |>
#   dplyr::select(seqnames,gene_id,gene_name) |>
#   left_join(dat,by=c(gene_id = "feature"))
# 
# # at first this looked great
# dat |>
#   filter(str_detect(label,"Sperm")) |>
#   ggplot(aes(MUT,WT)) +
#   ggdensity::geom_hdr_points() +
#   geom_abline(slope=1,intercept = 0,linetype="dashed",color="red") +
#   facet_wrap(~label,scales="free")
# 
# # doesnt seem to work very well because chr9 looks the same
# dat |>
#   filter(str_detect(seqnames,"chrX|chr9")) |>
#   ggplot(aes(seqnames,share.score)) +
#   geom_boxplot(outlier.shape = NA) +
#   #scale_y_sqrt() +
#   facet_wrap(~label)
