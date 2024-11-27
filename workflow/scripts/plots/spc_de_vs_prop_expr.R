Sys.setenv(R_PROFILE=".Rprofile")
source(Sys.getenv("R_PROFILE"))

library(tidyverse)
library(scater)
library(scuttle)
library(scran)

classifications <- read_tsv("data/dfam_classification.parsed.txt") |>
  filter(classification == "LINE")
coi <- "Spermatocyte"
de <- read_tsv("results/differential_expression/adult.tbl.pseudobulk.de.tsv.gz") |>
  filter(celltype == coi) |>
  mutate(class = case_when(FDR>0.05~"ns",
                           logFC > 0 ~"up",
                           logFC < 0 ~ "dn"))
#filter(feature %in% classifications$dfam_name)

sce <- read_rds("results/germ_cells/adult.sce.germ_cell.both_genotypes.subclustered.reintegrated.rds")

normcounts(sce) <- normalizeCounts(sce,
                log=F,
                size.factors = sizeFactors(sce),
                center_size_factors = F)

#sce <- sce[,sce$celltype == coi]
sce.wt <- sce[,sce$genotype =="WT"]
sce.mut <- sce[,sce$genotype =="MUT"]


entropy <- TSCAN::perCellEntropy(sce)

plotColData(sce,x="genotype",y=I(entropy),other_fields = "label") +
  facet_wrap(~label,nrow=2) +
  stat_summary(geom="point",color="red",fun = median) +
  ggpubr::stat_compare_means(size=3) +
  ylab("entropy")


df <- as_tibble(perFeatureQCMetrics(sce,assay.type="normcounts",subsets=list(wt=colnames(sce.wt),mut=colnames(sce.mut))),rownames="feature")

# "genes detected in fewer cells tend to be more differentially expressed"
df |>
  inner_join(de) |>
  ggplot(aes(detected,logFC,color=class)) +
  geom_point()


# "at a given expression level, de genes tend to be those not detected in as many (both genotype) cells"
df |>
  inner_join(de) |>
  arrange(logCPM) |>
  ggplot(aes(logCPM,detected,color=class)) +
  geom_point()

df |>
  inner_join(de) |>
  arrange(logCPM) |>
  ggplot(aes(logCPM,subsets_wt_detected,color=class)) +
  geom_point()

df |>
  inner_join(de) |>
  arrange(logCPM) |>
  ggplot(aes(logCPM,subsets_mut_detected,color=class)) +
  geom_point()


# 'genes that go up are expressed in more mutant cells than wt cells'
# lets say youre a low count gene because you're expressed in only a fraction of cells
# you need to be shared so all cells can use your products
# when bridges are gone, you concentrate in a few cells, thus reaching the detection
# threshold for inclusion in any given droplet's prepped cDNAs
# so perhaps losing intact bridges altered the dropout rate for these genes
df |>
  inner_join(de) |>
  dplyr::select(feature,subsets_wt_detected, subsets_mut_detected,FDR,logCPM,logFC,class) |>
  pivot_longer(-c(feature,FDR,logCPM,logFC,class)) |>
  mutate(bin=cut_interval(logCPM,length=1)) |>
  ggplot(aes(bin,value,fill=name)) +
  geom_boxplot() +
    facet_wrap(~class, labeller = label_both)


# next up - average count vs summed count vs pct detected vs significance
df |>
  inner_join(de) |>
  dplyr::select(feature,subsets_wt_detected, subsets_mut_detected,FDR,logCPM,logFC) |>
  ggplot(aes(subsets_wt_detected,subsets_mut_detected,color=factor(sign(logFC)),alpha=FDR<0.05)) +
  geom_point()



# expect that de genes are highly expressed in mut germs that that actually are expressing that gene
# but are detected in few wt cells and few mut cells and are lower expressed among expressing cells in wt 
mut_expr_mean <- normcounts(sce.mut) |>
  as.matrix() |>
  as_tibble(rownames="feature") |>
  pivot_longer(-feature) |>
  filter(value > 0) |>
  group_by(feature) |>
  summarize(subsets_expressed_mut_mean = mean(value))


wt_expr_mean <- normcounts(sce.wt) |>
  as.matrix() |>
  as_tibble(rownames="feature") |>
  pivot_longer(-feature) |>
  filter(value > 0) |>
  group_by(feature) |>
  summarize(subsets_expressed_wt_mean = mean(value))

expr_mean <- full_join(wt_expr_mean,mut_expr_mean) |>
  mutate(across(c(subsets_expressed_mut_mean,subsets_expressed_mut_mean),replace_na,0)) |>
  inner_join(de) |>
  inner_join(df)

# up genes have low mean expression, but high mean among expressing cells
expr_mean |>
  filter(class!="ns") |>
ggplot(aes((subsets_expressed_mut_mean + 1)/(subsets_mut_mean + 1), subsets_mut_detected,color=class)) +
  geom_point() +
  scale_y_log10()

expr_mean |>
  filter(class!="ns") |>
  ggplot(aes((subsets_expressed_mut_mean + 1)/(subsets_mut_mean + 1),(subsets_expressed_wt_mean + 1)/(subsets_wt_mean + 1),color=class)) +
  geom_point() +
  scale_y_log10() +
  scale_x_log10()

  
expr_mean |>
    #filter(class!="ns") |>
    ggplot(aes(class,(subsets_expressed_wt_mean + 1)/(subsets_wt_mean + 1),color=class)) +
    geom_boxplot() +
    scale_y_log10() +
    ggpubr::stat_compare_means()


expr_mean |>
  #filter(class!="ns") |>
  ggplot(aes(class,subsets_wt_detected,color=class)) +
  geom_boxplot() +
  ggpubr::stat_compare_means()






# gini



# expect that de genes are highly expressed in mut germs that that actually are expressing that gene
# but are detected in few wt cells and few mut cells and are lower expressed among expressing cells in wt 

library(DescTools)
mut_expr_gini <- logcounts(sce.mut) |>
  as.matrix() |>
  as_tibble(rownames="feature") |>
  filter(feature %in% de$feature) |>
  pivot_longer(-feature) |>
  group_by(feature) |>
  summarize(mut_gini = Gini(value))


wt_expr_gini <- logcounts(sce.wt) |>
  as.matrix() |>
  as_tibble(rownames="feature") |>
  filter(feature %in% de$feature) |>
  pivot_longer(-feature) |>
  group_by(feature) |>
  summarize(wt_gini = Gini(value))

TrajectoryUtils::

expr_gini <- full_join(wt_expr_gini,mut_expr_gini) |>
  inner_join(de) |>
  inner_join(df)


expr_gini



expr_gini |>
  filter(feature %in% classifications$dfam_name) |>
  ggplot(aes(wt_gini,mut_gini)) +
  geom_point() +
  geom_abline(intercept=0,slope=1)


expr_gini |>
  dplyr::select(feature,wt_gini,mut_gini,class) |>
  pivot_longer(-c(feature,"class")) |>
  ggplot(aes(name,value)) +
  geom_boxplot() +
  facet_wrap(~class) +
  ggpubr::stat_compare_means(paired=T)


expr_gini |>
  #filter(class!="ns") |>
  filter(feature %in% classifications$dfam_name) |>
  dplyr::select(feature,wt_gini,mut_gini,class,logFC) |>
  mutate(gini_ratio = log2(mut_gini/wt_gini)) |>
  ggplot(aes(gini_ratio,logFC,color=class)) +
  geom_point()
