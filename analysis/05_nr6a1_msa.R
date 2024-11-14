library(msa)
library(ggtree)
library(aaaView)
library(tidyverse)

args = commandArgs(trailingOnly=TRUE)

fasta_file <- args[1]
width <- args[2]
height <- args[3]
spacer <- as.numeric(args[4])
svg_out_file <- args[5]

input_aa <- Biostrings::readAAStringSet(fasta_file)
input_aa@ranges@NAMES <- gsub('nuclear receptor subfamily 6 group A member 1','NR6A1', input_aa@ranges@NAMES)
color_pick <- aaaView:::scheme_AA$Chemistry_AA
names(color_pick) <- row.names(aaaView:::scheme_AA)

msa_align <- msa(input_aa, method = 'ClustalOmega')

consensus <- AAStringSet(x=msaConsensusSequence(msa_align))
names(consensus) <- 'Consensus'
tidy_msa <- aaaView:::tidy_msa(c(msa_align@unmasked, consensus) )

# derive original coordinates for each protein
orig_position <- function(tidy_msa){
  data <- list()
  for (i in unique(tidy_msa$name)){
    data[[i]] <- tidy_msa %>% filter(name == i, character != '-') %>%
      mutate(orig_position = row_number())
  }
  tidy_msa %>% left_join(., data %>%
                           bind_rows() %>%
                           select(name, position, orig_position),
                         by = c('name','position')) %>%
    arrange(name, position)
}

tidy_msa <- orig_position(tidy_msa)
plot_data <- tidy_msa %>%
  mutate(group =
           substr(formatC(position, width = 5, format = "d", flag = "0"), 3, 3),
         Position =
           substr(formatC(position, width = 5, format = "d", flag = "0"), 4, 5)) %>%
  mutate(Protein = case_when(name != 'Consensus' ~ name,
                             TRUE ~ 'Consensus'),
         Protein = factor(Protein, levels = rev(c("Consensus", input_aa@ranges@NAMES)))) %>%
  dplyr::rename(AA = character)
plot_data$group <- factor(plot_data$group)
levels(plot_data$group) = paste0(plot_data$group, " (Hundreds Position)") %>% unique()


svg(svg_out_file, width = as.integer(width), height = as.integer(height))
main_plot <- plot_data %>%
  ggplot(aes(x=Position,y=Protein, label = AA, fill = AA)) +
  facet_wrap(~group,  ncol = 1) +
  geom_tile() +
  geom_text(size =3) +
  theme_minimal() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())  +
  scale_fill_manual(values = color_pick) +
  xlab("") + 
  ylab("") + 
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    axis.text.x = element_blank(),
    strip.text.y = element_blank()) + 
  theme(legend.position="none") +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

facet_end_coords <- bind_rows(plot_data %>% 
                                filter(position %in% seq(0,1e6,100)) %>% 
                                mutate(group = paste0(
                                  (str_extract(group, "\\d") %>% as.integer())-1,
                                  " (Hundreds Position)")),
                              plot_data %>% as_tibble() %>% 
                                group_by(name) %>% slice_max(orig_position))%>% 
  ggplot(aes(x=0,y=Protein, label = orig_position)) +
  facet_wrap(~group,  ncol = 1) +
  geom_text(size = 3) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())  +
  scale_fill_manual(values = color_pick) +
  xlab("") + 
  ylab("") + 
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    strip.text.y = element_blank()) +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
cowplot::plot_grid(main_plot , NULL, nrow = 1,
                   facet_end_coords, rel_widths = c(1,spacer, 0.3), align = 'hv', axis = "bt") 
dev.off()


# 
# p <- ggtree(tree_seq_nwk) + geom_tiplab(size=3)
# msaplot(p, AA_sequence, offset=3, width=2)
