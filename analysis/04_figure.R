# 2025 01 28
# update for resubmission

# 2025 03 27
# update final (?) submission
# remove Aldiri dataset from fetal plot
# and add sample counts to boxplot
library(tidyverse)
load('data/01_image.Rdata')
load('data/02_fetal_plot.Rdata')
load('data/02_eyeIntegration.files.Rdata')
load('data/03_correlation.image.Rdata')
body_plot <- body_data %>% 
  ggplot(aes(x=Tissue, y=log1p(CPM))) + 
  geom_boxplot() + scale_x_discrete(guide = guide_axis(angle = 90)) +
  geom_text(data = .  %>% group_by(Tissue) %>% summarise(Count = n()),
            aes(label = Count, y = 7)) +
  cowplot::theme_cowplot() + coord_cartesian(ylim=c(0,8)) + 
  ylab("log1p(bulk\nNR6A1)") + xlab("GTEx Bulk Tissue")

eiad_plot <- eiad_data %>% 
  mutate(Tissue = case_when(Tissue == 'ESC' ~ 'ESC, iPSC',
                            Tissue == 'iPSC' ~ 'ESC, iPSC',
                            Tissue == 'Trabecular Meshwork' ~ 'Trabecular\nMeshwork',
                            TRUE ~ Tissue)) %>% 
  ggplot() + 
  geom_boxplot(aes(x=Tissue, y=log1p(CPM))) +  
  ylab("log1p(bulk\nNR6A1)") +
  geom_text(data = .  %>% group_by(Tissue) %>% summarise(Count = n()) %>% 
              mutate(CPM = c(6, 8, 9.5, 8, 8, 8, 9.5, 9.5)),
            aes(x=Tissue, y=(CPM), label = Count)) +
  xlab("eyeIntegration Bulk Tissues") +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  cowplot::theme_cowplot() + 
  coord_cartesian(ylim=c(0,9.5))


hrca <- pb_long %>% 
  left_join(conv_table %>% group_by(ENSEMBL) %>% summarise(SYMBOL = paste(SYMBOL, collapse= ', ')), 
            by = c("name" = "ENSEMBL")) %>% 
  filter(SYMBOL == 'NR6A1') %>% 
  ggplot(aes(x=class,y=(value))) + 
  geom_boxplot(outlier.colour = NA) +
  ggbeeswarm::geom_quasirandom(size = 0.5) +
  geom_text(aes(label = Count),data = .  %>% group_by(class) %>% summarise(Count = n()) %>% 
              mutate(value = rep(9,10))) +
  geom_text(aes(label = Mean),data = .  %>% group_by(class) %>% summarise(Mean = mean(value) %>% 
                                                                            round(., digits = 1)) %>% 
              mutate(value = rep(8.5,10))) +
  geom_text(aes(label = Median),data = .  %>% group_by(class) %>% summarise(Median = median(value) %>% 
                                                                              round(., digits = 1)) %>% 
              mutate(value = rep(8,10))) +
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
                            direction = 'y', max.overlaps = Inf, force = 100) +
  xlab("Correlation of gene \nexpression with NR6A1 across\nFetal Retina and RPE") +
  cowplot::theme_cowplot() + ylab('') 

eiad_eye_fetal %>% 
  mutate(Coloboma = case_when(Gene %in% colo$Gene ~ 'Coloboma')) %>% 
  select(-Gene2N) %>% 
  write_csv("nr6a1_correlation_table.csv")
top <- cowplot::plot_grid(eiad_plot, fetal_plot, cor + ylab("density"), nrow =1, 
                          labels = 'auto', rel_widths = c(0.9,1.2,1.1))

# svg('nr6a_expression_figure.20250128.svg', width = 10, height = 7)
# cowplot::plot_grid(top, 
#                    NULL,
#                    body_plot, 
#                    nrow = 3, rel_heights = c(1,-0.35,0.8), 
#                    align = 'h',
#                    #rel_widths = c(1,1,1),
#                    labels = c('','','d'), label_size = 12)
# dev.off()


pdf('nr6a_expression_figure.20250402.01.pdf', width = 16, height = 9)
cowplot::plot_grid(top,
                   NULL,
                   body_plot,
                   nrow = 3, rel_heights = c(0.9,-0.3,0.8),
                   align = 'h',
                   #rel_widths = c(1,1,1),
                   labels = c('','','d'), label_size = 11.5)
dev.off()

# data table out



# pdf("nr6a_single_cell_supplemental_expression_figure.20250402.pdf", width = 4, height = 3)
# hrca
# dev.off()
