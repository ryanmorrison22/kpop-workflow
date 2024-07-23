#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggtree))

set.seed(1)

args = commandArgs(trailingOnly=TRUE)
if (length(args)!=3) {
  stop("Three arguments must be supplied - <original_nwk> <updated_nwk> <output_file", call.=FALSE)
}

# Set arguments as variables
original_nwk_file = args[1]
updated_nwk_file = args[2]
output_file = args[3]

original_tree_file<-read.tree(original_nwk_file)
updated_tree_file<-read.tree(updated_nwk_file)

samples2highlight <- updated_tree_file$tip.label[!(updated_tree_file$tip.label %in% original_tree_file$tip.label)]

colours = as.data.frame(updated_tree_file$tip.label) %>% mutate(group = ifelse(!(updated_tree_file$tip.label %in% original_tree_file$tip.label), "updated", "original"))

original_tree <- ggtree(original_tree_file, layout="circular",color="grey75") +
  geom_tiplab() +
  scale_x_continuous(expand=expansion(0.5)) + 
  ggtitle("Original Tree") +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))

updated_tree <- ggtree(updated_tree_file, layout="circular",color="grey75") %<+% colours +
  geom_tiplab(aes(colour = factor(group)), show.legend =F) + 
  scale_colour_manual(values = c("black", "red"), na.translate=FALSE, name = "Group") +
  scale_x_continuous(expand=expansion(0.5)) + 
  ggtitle("Updated Tree") +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))
  
tree <- original_tree +
  updated_tree 

ggsave(file = output_file, plot = tree, height = 20, width = 20)
