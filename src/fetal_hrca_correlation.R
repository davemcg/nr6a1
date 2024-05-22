library(tidyverse)
library(matrixStats)
library(WGCNA)
library(SummarizedExperiment)
library(spqn)

pb <- data.table::fread("data/chen_fetal_hrca.pb_raw.csv.gz")

pb_mat <- pb[,2:ncol(pb)] %>% as.matrix()
row.names(pb_mat) <- pb %>% pull(1)
pb_mat <- pb_mat[pb_mat %>% rowSums() %>% log1p() > 12,]


pb_transform <- metamoRph::normalize_data(t(pb_mat), log1p = TRUE)
pb_transform_g <-  pb_transform[rowMedians(pb_transform) > 0,]
                 
ave_logrpkm <- rowMeans(pb_transform_g)
logrpkm <- pb_transform_g - ave_logrpkm # mean centering
logrpkm  <- logrpkm / matrixStats::rowSds(logrpkm) # variance scaling

# remove PCs from the gene expression matrix after scaling each gene to have mean=0 and variance=1
logrpkm_4pc <- removePrincipalComponents(t(scale(t(logrpkm))), n = 4)

# base correlation
cor_m <- cor(t(logrpkm_4pc))

# more correlations with higher expression
plot_signal_condition_exp(cor_m, ave_logrpkm, signal=0)
# do norm
cor_m_spqn <- normalize_correlation(cor_m, ave_exp=ave_logrpkm, ngrp=20, size_grp=300, ref_grp=18)
row.names(cor_m_spqn) <- row.names(cor_m)
colnames(cor_m_spqn) <- colnames(cor_m)
# balanced correlation results
plot_signal_condition_exp(cor_m_spqn, ave_logrpkm, signal=0)

# make it long
cor_long <- cor_m_spqn %>% as_tibble(rownames = 'ensembl1') %>% 
  pivot_longer(-ensembl1, names_to = 'ensembl2', values_to = 'correlation')

# 
conv_table <- AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, 
                                    keys=colnames(pb),
                                    columns=c("ENSEMBL","SYMBOL", "MAP","GENENAME", "ENTREZID"), keytype="ENSEMBL") 
# annotate with human readable gene names
# remove low correlations
nr6a1_cor <-  cor_long %>% 
  filter(ensembl1 == 'ENSG00000148200') %>% 
  filter(ensembl1 != ensembl2,
         ensembl2 != 'NA') %>% 
  left_join(conv_table %>% group_by(ENSEMBL) %>% 
              summarise(Gene2 = paste0(SYMBOL, collapse = ', ')),
            by = c("ensembl2" = "ENSEMBL")) %>% 
  arrange(-abs(correlation)) %>% 
  mutate(rank = row_number())
  
# rank of known coloboma genes
colo <- read_csv("https://raw.githubusercontent.com/davemcg/eyeMarkers/master/lists/george_brooks_coloboma_2020.csv")
nr6a1_cor %>% 
  mutate(Coloboma = case_when(Gene2 %in% colo$Gene ~ 'Coloboma')) %>% 
  filter(!is.na(Coloboma))

write_tsv(nr6a1_cor, file = 'data/hrca_fetal_nr6a1_cor.tsv.gz')
