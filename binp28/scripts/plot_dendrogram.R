#!/usr/bin/env Rscript

rm(list=ls())

# Set CRAN and Bioconductor repositories
options(repos = c(CRAN = "https://cloud.r-project.org/"))

# Install Bioconductor manager if not installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# List required packages
packages <- c("ape", "ggplot2", "tidyverse", "ggtree")

# Install missing packages
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    if (pkg == "ggtree") {
      BiocManager::install(pkg)
    } else {
      install.packages(pkg)
    }
    library(pkg, character.only = TRUE)
  }
}

library(ape)
library(ggtree)

print("All packages loaded successfully!")

# Parse command-line arguments
args = commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript phylo_plot.R <tree_file> <outgroup> <img_path>")
}

tree_file = args[1]
outgroup = args[2]
save_img = args[3]

# Check if tree file exists
if (!file.exists(tree_file)) {
  stop(paste("Tree file not found:", tree_file))
}

# Read and root the tree
tree <- read.tree(tree_file)
if (!(outgroup %in% tree$tip.label)) {
  stop(paste("Outgroup", outgroup, "not found in the tree!"))
}
rooted_tree <- root(tree, outgroup = outgroup, resolve.root = TRUE)
rooted_tree <- ladderize(rooted_tree)  # Improve readability

# Save the tree plot
png(save_img, width = 1600, height = 1000, res = 200)

ggtree(rooted_tree, layout = "rectangular") + 
  geom_tiplab(size = 3, hjust = -0.1) + 
  theme_tree2() +  
  xlim(0, max(node.depth.edgelength(rooted_tree)) * 1.2) +  
  labs(title = "Rooted Phylogenetic Tree")

dev.off()

print("Done!")


