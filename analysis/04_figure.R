# 2025 01 28
# update for resubmission
library(tidyverse)
load('data/01_image.Rdata')
load('data/02_fetal_plot.Rdata')
load('data/02_eyeIntegration.files.Rdata')
load('data/03_correlation.image.Rdata')
body_plot <- body_data %>% 
  ggplot(aes(x=Tissue, y=log1p(CPM))) + 
  geom_boxplot() + scale_x_discrete(guide = guide_axis(angle = 90)) +
  cowplot::theme_cowplot() + coord_cartesian(ylim=c(0,8)) + 
  ylab("log1p(bulk\nNR6A1)") + xlab("GTEx Bulk Tissue")

eiad_plot <- eiad_data %>% 
  mutate(Tissue = case_when(Tissue == 'ESC' ~ 'ESC, iPSC',
                            Tissue == 'iPSC' ~ 'ESC, iPSC',
                            Tissue == 'Trabecular Meshwork' ~ 'Trabecular\nMeshwork',
                            TRUE ~ Tissue)) %>% 
  ggplot(aes(x=Tissue, y=log1p(CPM))) + 
  geom_boxplot() +  ylab("log1p(bulk NR6A1)") + 
  xlab("eyeIntegration Bulk Tissues") +
scale_x_discrete(guide = guide_axis(angle = 90)) +
  cowplot::theme_cowplot() + 
  coord_cartesian(ylim=c(0,8))


hrca <- pb_long %>% 
  left_join(conv_table %>% group_by(ENSEMBL) %>% summarise(SYMBOL = paste(SYMBOL, collapse= ', ')), 
            by = c("name" = "ENSEMBL")) %>% 
  filter(SYMBOL == 'NR6A1') %>% 
  ggplot(aes(x=class,y=(value))) + 
  geom_boxplot(outlier.colour = NA) +
  ggbeeswarm::geom_quasirandom(size = 0.5) +
  ylab("log1p(snRNA\nNR6A1)") + xlab("HRCA Pseudobulk Cell Types") +
  scale_x_discrete(guide = guide_axis(angle = 90)) + 
  cowplot::theme_cowplot() 

cor <- eiad_eye_fetal %>% 
  mutate(Coloboma = case_when(Gene %in% colo$Gene ~ 'Coloboma'))  %>% 
  ggplot(aes(x=correlation)) + 
  geom_density() + 
  geom_point(aes(x=correlation,y=0), 
             data = . %>% filter(!is.na(Coloboma)), size = 1) +
  ggrepel::geom_label_repel(aes(x=correlation, y=0, label = Gene),
                            data = . %>% filter(!is.na(Coloboma)) %>% head(10),
                            direction = 'y', max.overlaps = Inf) +
  xlab("Correlation of NR6A1\nwith eyeIntegration Fetal\nRetina and RPE") +
  cowplot::theme_cowplot() + ylab('') 

top <- cowplot::plot_grid(fetal_plot, hrca, cor, nrow =1, labels = 'auto', rel_widths = c(1.3,0.9,1))

svg('nr6a_expression_figure.20250128.svg', width = 10, height = 7)
cowplot::plot_grid(top, 
                   NULL,
                   body_plot, 
                   nrow = 3, rel_heights = c(1,-0.35,0.8), 
                   align = 'h',
                   #rel_widths = c(1,1,1),
                   labels = c('','','d'), label_size = 12)
dev.off()

