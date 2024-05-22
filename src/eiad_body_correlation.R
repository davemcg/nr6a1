args = commandArgs(trailingOnly=TRUE)

library(tidyverse)
library(matrixStats)
library(WGCNA)
library(SummarizedExperiment)
library(spqn)

gene_counts <- data.table::fread("~/git/EiaD_build/salmon_counts/gene_counts.csv.gz")
metadata <- data.table::fread("~/git/EiaD_build/data/eyeIntegration23_meta_2023_09_01.built.csv.gz")
gene_counts_mat <- gene_counts[,2:ncol(gene_counts)] %>% as.matrix()
row.names(gene_counts_mat) <- gene_counts %>% pull(1)

metadata <- metadata %>% filter(sample_accession %in% colnames(gene_counts))

i <- args[1]
#for (i in (metadata %>% filter(Cohort == 'Body') %>% pull(Tissue) %>% unique())[16:31]){
print(i)
set.seed(1243)

samps <- (metadata %>% 
            filter(Cohort == 'Body') %>% 
            ungroup() %>% 
            filter(Tissue == i) %>% 
            sample_n(50, replace = TRUE) %>% 
            pull(sample_accession)) %>% 
  unique()

samps <- samps[samps %in% colnames(gene_counts_mat)]
gene_counts_mat_loop <- gene_counts_mat[,samps]
print(dim(gene_counts_mat_loop))
gc_transform <- metamoRph::normalize_data((gene_counts_mat_loop), log1p = TRUE)
gc_transform_g <-  gc_transform[rowMedians(gc_transform) > 0,]

ave_logrpkm <- rowMeans(gc_transform_g)
logrpkm <- gc_transform_g - ave_logrpkm # mean centering
logrpkm  <- logrpkm / matrixStats::rowSds(logrpkm) # variance scaling

# remove PCs from the gene expression matrix after scaling each gene to have mean=0 and variance=1
logrpkm_4pc <- removePrincipalComponents(t(scale(t(logrpkm))), n = 4)

# base correlation
cor_m <- cor(t(logrpkm_4pc))

# more correlations with higher expression
#plot_signal_condition_exp(cor_m, ave_logrpkm, signal=0)
# do norm
cor_m_spqn <- normalize_correlation(cor_m, ave_exp=ave_logrpkm, ngrp=20, size_grp=300, ref_grp=18)
row.names(cor_m_spqn) <- row.names(cor_m)
colnames(cor_m_spqn) <- colnames(cor_m)
# balanced correlation results
#plot_signal_condition_exp(cor_m_spqn, ave_logrpkm, signal=0)

# make it long
cor_long <- cor_m_spqn %>% as_tibble(rownames = 'Gene1') %>% 
  pivot_longer(-Gene1, names_to = 'Gene2', values_to = 'correlation')

# nr6a1 specific
nr6a1_cor <- cor_long %>% 
  filter(grepl('NR6A1', Gene1), Gene1 != Gene2) %>% 
  arrange(-abs(correlation)) %>% 
  mutate(rank = row_number()) %>% 
  mutate(Gene2N = gsub(" \\(.*","",Gene2)) %>% 
  mutate(Tissue = i)

write_tsv(nr6a1_cor, file = paste0('data/eyeIntegration_', i, '_nr6a1_cor.tsv.gz'))
#}
