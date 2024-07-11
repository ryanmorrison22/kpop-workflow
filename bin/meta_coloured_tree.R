#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggtree))
suppressPackageStartupMessages(library(ggtreeExtra))
suppressPackageStartupMessages(library(ggnewscale))
suppressPackageStartupMessages(library(RColorBrewer))

set.seed(1)

args = commandArgs(trailingOnly=TRUE)
if (length(args)!=3) {
  stop("Three arguments must be supplied - <meta data file> <nwk file> <output file name>", call.=FALSE)
}

# Set arguments as variables
phen_file = args[1]
nwk_file = args[2]
output_fileName = args[3]

# Load file contents and remove file name extensions
tree_file<-read.tree(nwk_file) 
phen<-read.table(phen_file,
                 header=TRUE,sep="\t", fill=TRUE)
phen$fileName <- sub('\\.fasta$|\\.fasta.gz$|\\.fastq$|\\.fastq.gz$', '', phen$fileName)

# Set the colours to be linked to the levels so they match across plots
custom_palette = c("#e40046", "#84bd00", "#ffb81c", "#00a5df", "#ff7f32", "#f032e6", "#fabed4", "#dcbeff", "#aaffc3", "#582c83", "#a9a9a9", "#007c91", "red", "#9A6324", "#d5cb9f", "#454545") 
colour_palettes = c("Accent", "Dark2", "Paired", "Set1", "Set2", "Set3", "Pastel1", "Pastel2", custom_palette) 

# Generate the tree with leaves colored by meta data
if ("class" %in% colnames(phen)) {
  colours = colorRampPalette(custom_palette)(n_distinct(phen$class))
  col = setNames(colours, str_sort(levels(factor(phen$class)), numeric = TRUE))
  tree <- ggtree(tree_file, layout="circular",color="grey75") %<+% phen + 
    geom_tippoint(pch=16, size = 5, aes(colour = factor(class, levels = str_sort(levels(factor(phen$class)), numeric = TRUE)))) +
    scale_colour_manual(values = col, na.translate=FALSE) +
    guides(color = guide_legend(title = "Class")) +
    theme(legend.title = element_text(size=40), legend.text = element_text(size=30))
} else {
  tree <- ggtree(tree_file, layout="circular",color="grey75") %<+% phen + 
    geom_tiplab() +
    theme(legend.title = element_text(size=40), legend.text = element_text(size=30))
}

extra_columns = phen %>% select(-any_of(c("fileName", "class")))

if (ncol(extra_columns) > 1) {
  count = 1
  for (i in colnames(extra_columns)){
    tree <- tree + 
      new_scale_fill() +
      geom_fruit(geom=geom_tile,
        mapping=aes(fill=factor(.data[[i]], levels = str_sort(levels(factor(.data[[i]]))))),      
        offset = 0.25,width = 0.25) + 
      scale_fill_brewer(palette=colour_palettes[count]) +
      labs(fill = i)
    count = count + 1
    if (count > 9) { #reset the count to stay within palette range
      count = 1
    }
  }
}

# Save plot
ggsave(file = output_fileName, plot = tree, height = 40, width = 40)
