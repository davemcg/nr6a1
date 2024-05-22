library(tidyverse)
library(ggtext)
metadata <- data.table::fread("~/git/EiaD_build/data/eyeIntegration23_meta_2023_09_01.built.csv.gz")
gene_counts <- data.table::fread("~/git/EiaD_build/salmon_counts/gene_counts.csv.gz")
gene_counts_mat <- gene_counts[,2:ncol(gene_counts)]
gene_cpm_mat <- edgeR::cpm(gene_counts_mat)
row.names(gene_cpm_mat) <- gene_counts$Gene

gene_counts_long <- gene_cpm_mat %>%
  as_tibble(rownames = 'Gene') %>%
  pivot_longer(cols=2:ncol(gene_counts), names_to = 'sample_accession', values_to = 'CPM')
gene_counts_long <- left_join(gene_counts_long,
                              metadata %>%
                                select(sample_accession, Tissue, Sub_Tissue, Perturbation,
                                       Source, Age, Sex_ML, study_title) %>%
                                unique(),
                              by = 'sample_accession')
plot_data <- gene_counts_long %>%
  #filter(Tissue %in% c("Cornea","RPE","Retina","Liver","Kidney")) %>%
  mutate(ID = paste(Tissue, Source, Age, sep = ' - ')) %>%
  filter(!grepl("AMD","Perturbation")) %>% 
  mutate(gene_id = gsub(" \\(.*","",Gene)) %>% 
  filter(gene_id %in% "NR6A1")

body_data <- plot_data %>% filter(grepl("GTEx",study_title))
eiad_data <- plot_data %>% filter(!grepl("GTEx",study_title),
                                  !Tissue %in% c("EyeLid"),
                                  !grepl("AMD", Perturbation))

plot <- plot_data %>%
  mutate(Sub_Tissue = case_when(is.na(Sub_Tissue) ~ '', TRUE ~ Sub_Tissue),
         Source = case_when(is.na(Source) ~ '', TRUE ~ Source),
         Age = case_when(is.na(Age) ~ '', TRUE ~ Age)) %>%
  mutate(Sub_Tissue = glue::glue("<span style='color:#000000FF'>{Sub_Tissue}</span>"),
         Source = glue::glue("<span style='color:#1E46A2FF'>{Source}</span>"),
         Age = glue::glue("<span style='color:#FB323BFF'>{Age}</span>")) %>%
  mutate(Sex_ML = case_when(is.na(Sex_ML) ~ "Unknown",
                            TRUE ~ Sex_ML),
         ID = gene_id) %>%
  ggplot(data=.,aes(x=interaction(Source, Sub_Tissue, Age, sep = ' | '),y=log2(CPM+1),
                    color = Tissue,
                    fill = Tissue)) +
  #geom_violin(alpha=0.5, scale = 'width') +
  geom_boxplot(alpha=0.7, outlier.shape = NA, width = 0.6, fill = 'black') +
  cowplot::theme_cowplot(font_size = 15) + theme(axis.text.x = element_text(angle = 90, hjust=1, vjust = 0.2)) +
  ggtitle('Pan-Human NR6A1 Expression') +
  ylab("log2(CPM+1)") +
  scale_shape_manual(values=c(0:2,5,6,15:50)) +
  theme(strip.background = element_rect(fill = 'black'),
        strip.text = element_text(color = 'white'),
        panel.background = element_rect(fill = 'gray90'),
        plot.margin=grid::unit(c(0,0,0,0.1), "cm"),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.key.size= unit(0.2, "cm"),
        legend.spacing = unit(0.2, "cm")) +
  coord_flip() +
  facet_grid(rows = vars(Tissue), cols = vars(ID), scales = 'free', space = 'free') +
  theme(strip.text.y.right = element_text(angle = 0)) +
  theme(
    axis.text.y = element_markdown(),
    axis.title.y = element_markdown()) +
  labs(x = "<span style='color:#1E46A2FF'>Source</span> |
<span style='color:#000000FF'>Sub Tissue</span> |
<span style='color:#FB323BFF'>Age</span>")

svg(filename = '~/git/nr6a1/eyeintegration_NR6A1.svg', height = 20, width = 10)
plot
dev.off()

save(plot_data, body_data, eiad_data, plot, file = 'data/02_eyeIntegration.files.Rdata')
