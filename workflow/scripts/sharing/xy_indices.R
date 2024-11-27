Sys.setenv(R_PROFILE=".Rprofile")
source(Sys.getenv("R_PROFILE"))

library(tidyverse)
library(scater)
library(scran)
library(scuttle)

sce <- read_rds(ifelse(exists("snakemake"),snakemake@input$sce,
                       "results/germ_cells/adult.sce.germ_cell.both_genotypes.subclustered.reintegrated.rds"))


# we can detect initiation of msci in spc, and consequences of bridge loss in spd
colData(sce) |>
  as_tibble(rownames="cell") |>
  dplyr::select(cell,genotype,label,total,nexprs,matches("chr") & contains("detected")) |>
  pivot_longer(-c(cell,genotype,label,total,nexprs)) |>
  filter(str_detect(name,"chrX|chrY|chr9_")) |>
  ggplot(aes(label,value+1,fill=name)) +
  #geom_boxplot(outlier.shape = NA) +
  geom_violin(scale = "width",adjust=1) +
  scale_y_log10() +
  facet_wrap(genotype~name,scales = "free_y") +
  theme(axis.text.x = element_text(angle=45,hjust=1))

# another look that makes mutual exclusivity clear
plotColData(sce,x = "subsets_chrX_detected",y = "subsets_chrY_detected",other_fields = c("celltype","genotype","label")) +
  facet_wrap(genotype~label,nrow=2,scales="free")

# chr9 as control
plotColData(sce,x = "subsets_chrX_detected",y = "subsets_chr9_detected",other_fields = c("celltype","genotype","label")) +
  facet_wrap(genotype~label,nrow=2,scales="free")


# a defect in msci?
# "persists in haploid spermatids"-https://www.nature.com/articles/s41467-022-34295-5
colData(sce) |>
  as_tibble(rownames="cell") |>
  dplyr::select(cell,genotype,label,total,nexprs,matches("chr") & contains("detected")) |>
  pivot_longer(-c(cell,genotype,label,total,nexprs)) |>
  filter(str_detect(name,"chrX|chrY|chr9_")) |>
  ggplot(aes(label,value,color=genotype)) +
  geom_boxplot(outlier.shape = NA) +
  #geom_point(position = position_jitterdodge(jitter.width = 0.2)) +
  #geom_violin(scale = "width") +
  facet_wrap(~name,scales = "free_y") +
  theme(axis.text.x = element_text(angle=45,hjust=1)) +
  scale_color_grey(end = .6) +
  scale_y_log10()


# te expressers are predominantly ones not expressing the Y
# L1s are enriched on the X
# could be a factor on the y required for msci?
# zfy1/2 are required for msci ~maintenance~ because it initiates when cells are
# not haploid
colData(sce) |>
  as_tibble(rownames="cell") |>
  dplyr::select(cell,genotype,label,subsets_TE_percent,subsets_TE_sum,total,nexprs,matches("chr")) |>
  group_by(genotype,label) |>
  #pivot_longer(-c(cell,genotype,label,total,nexprs,subsets_TE_percent)) |>
  #group_by(genotype,label,name) |>
  #arrange(genotype,label,name,-value) |>
  #mutate(value = row_number()) |>
  #ungroup() |>
  #pivot_wider(names_from = name,values_from = value) |>
  filter(str_detect(label,"Spermatid|Spermatocyte")) |>
  arrange(subsets_TE_percent) |>
  ggplot(aes(subsets_chrX_detected,subsets_chrY_detected, color=scales::squish(subsets_TE_percent,range=c(0,15)))) +
  geom_point() +
  #ggdensity::geom_hdr_points(size=0.2) +
  facet_wrap(genotype~label,nrow=2,scales="free") +
  theme_dark() +
  scale_color_viridis_c()



# need fig showing total umi count similar for 
# x and y-bearing mutant spermatids


y.linked <- rowRanges(sce) |>
  as_tibble() |>
  filter(seqnames=="chrY") |>
  dplyr::select(gene_name,gene_id) |>
  filter(gene_id %in% rownames(sce)[rowSums(counts(sce)) > 20])

read_tsv("results/differential_expression/adult.tbl.germ_cell.reprocessed.de.tsv.gz") |>
  filter(FDR < 0.05 & feature %in% y.linked$gene_id) |>
  left_join(y.linked,by=c(feature="gene_id")) |>
  print(n=Inf)


