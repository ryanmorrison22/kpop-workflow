#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("tidyverse"))

set.seed(1)

args = commandArgs(trailingOnly=TRUE)
if (length(args)!=4) {
  stop("Four arguments must be supplied - <embedding file> <num of clusters> <centroid/kmeans> <output file name>", call.=FALSE)
}

# Set arguments as variables
embed_file = args[1]
cluster_num = args[2]
clust_type = args[3]
output_fileName = args[4]

embeddings <- read.table(embed_file, header = TRUE, row.names = 1, sep = '\t')
dis_mat <- dist(embeddings, method = "euclidean")

if (clust_type == "centroid") {
  ## Centroid
  centroid_hclust <- hclust(dis_mat, method = "centroid")
  i <- as.numeric(cluster_num)
  while (i > 1) { # If the cluster number is less than the number of dimensions then choose highest number of valid clusters
    skip_to_next <- FALSE
    tryCatch(c_clusters <- as.data.frame(cutree(centroid_hclust, k = i)), error = function(e) { skip_to_next <<- TRUE})
    i <- i - 1
    if(skip_to_next) { next } else { break }
  }
  c_clusters <- rownames_to_column(c_clusters, "fileName")
  colnames(c_clusters) <- c("fileName", "class")
  write.table(file = output_fileName, c_clusters, quote = F, row.names = F, sep = '\t')
}

if (clust_type == "kmeans") {
  ## K-means
  i <- as.numeric(cluster_num)
  while (i > 1) { # If the cluster number is less than the number of dimensions then choose highest number of valid clusters
    skip_to_next <- FALSE
    tryCatch(kmeans_clust <- kmeans(dis_mat, centers = i), error = function(e) { skip_to_next <<- TRUE})
    i <- i - 1
    if(skip_to_next) { next } else { break }
  }
  k_clusters <- as.data.frame(kmeans_clust$cluster)
  k_clusters <- rownames_to_column(k_clusters, "fileName")
  colnames(k_clusters) <- c("fileName", "class")
  write.table(file = output_fileName, k_clusters, quote = F, row.names = F, sep = '\t')
}
