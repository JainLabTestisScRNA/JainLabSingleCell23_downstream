library(tidyverse)
library(patchwork)
library(clusterProfiler)

# ------------------------------------------------------------------------------
g_k_choice <- read_rds("results/tradeseq/juvenile_13d_wt_null.cellranger.germ_cell.g_k_choice.rds")

# ------------------------------------------------------------------------------
# smoothed timecourse + clustering
res <- read_rds("results/tradeseq/juvenile_13d_wt_null.cellranger.germ_cell.smoothed_cluster_timecourse.rds")

ord <- res |> filter(condition=="WT") |>
  group_by(cluster, pseudotime.bin) |>
  summarise(yhatScaled = mean(yhatScaled)) |>
  group_by(cluster) |>
  slice_max(yhatScaled, n=1) |>
  arrange(pseudotime.bin) |>
  pull(cluster)

res <- res |> mutate(cluster = fct_relevel(cluster, ord))

theme_set(theme_bw())

g_smooth <- res |>
  ggplot(aes(pseudotime.bin,gene_name,fill=yhatScaled)) +
  geom_tile() +
  scale_fill_viridis_c(begin = 0.2, end= 1, option = "C",) +
  facet_grid(cluster~condition, scales="free",space = "free") +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  scale_x_continuous(expand = expansion())

g_diff <- res |>
  filter(condition == "wt") |>
  mutate(condition = "lfcMutant.WT") |>
  mutate(gene_name = fct_reorder(gene_name, lfcMutant.WT)) |>
  ggplot(aes(pseudotime.bin,gene_name,fill= lfcMutant.WT)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red") +
  facet_grid(cluster~., scales="free",space = "free") +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

g_pattern <- res |>
  #filter(condition=="wt") |>
  group_by(condition, cluster,pseudotime.bin) |>
  summarise(mean.yhatScaled = mean(yhatScaled), .groups = "drop") |>
  ggplot(aes(pseudotime.bin,mean.yhatScaled, color=condition)) +
  geom_path() +
  facet_grid(cluster~.)

g_smooth + g_diff + g_pattern + plot_layout(widths = c(2,1,0.5)) & ylab("gene") & xlab("pseudotime bin") & scale_x_continuous(expand = expansion()) & theme(legend.position = "bottom")

theme_set(theme_classic())

g_n_tes <- res |>
  filter(!str_detect(feature,"ENSMUSG")) |>
  dplyr::select(feature, cluster) |>
  distinct()  |>
  ggplot(aes(cluster)) +
  geom_bar() + ylab("n TEs")

g_n_l1 <- res |>
  filter(str_detect(feature,"^L1")) |>
  dplyr::select(feature, cluster) |>
  distinct()  |>
  ggplot(aes(cluster)) +
  geom_bar() + ylab("n LINE-1s")

g_n_tes / g_n_l1

# ------------------------------------------------------------------------------
# clusterprofiler
library(clusterProfiler)
clp <- read_rds("results/tradeseq/juvenile_13d_wt_null.cellranger.germ_cell.clusterprofiler.rds")

clp$per.cluster |> dotplot()

clp$per.cluster.direction |> dotplot()

clp$per.cluster.direction |>
  filter(p.adjust < 0.01) |>
  dotplot()
