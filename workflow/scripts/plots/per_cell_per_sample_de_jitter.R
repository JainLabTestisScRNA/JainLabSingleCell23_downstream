library(tidyverse)
library(scater)
library(scran)
library(scuttle)

sce <- read_rds("results/integration/adult.sce.integrated.clustered.celltypes.rds")

de <- read_tsv("results/differential_expression/adult.tbl.pseudobulk.de.tsv.gz")

Biostrings::readDNAStringSet("data/dfam_sequences.fasta")

filter(de,!str_detect(feature,"ENSMUSG"))
