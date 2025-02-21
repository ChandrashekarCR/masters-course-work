#!/usr/bin/env Rscript

rm(list = ls())

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
  stop("Usage: Rscript phylo_plot_bootstrap.R <tree_file> <outgroup> <img_path>")
}

tree_file = args[1]
outgroup = args[2]
save_img = args[3]

# Check if tree file exists
if (!file.exists(tree_file)) {
  stop(paste("Tree file not found:", tree_file))
}

# Read and root the tree (handling potential polytomies)
tree <- read.tree(tree_file)

if (!(outgroup %in% tree$tip.label)) {
  stop(paste("Outgroup", outgroup, "not found in the tree!"))
}

rooted_tree <- root(tree, outgroup = outgroup, resolve.root = TRUE)
rooted_tree <- ladderize(rooted_tree) # Improve readability



# Create the plot
png(save_img, width = 1600, height = 1000, res = 200)

p <- ggtree(rooted_tree, layout = "rectangular") +
  geom_tiplab(size = 3, hjust = -0.1) +
  theme_tree2() +
  xlim(0, max(node.depth.edgelength(rooted_tree)) * 1.2) +
  labs(title = "Rooted Phylogenetic Tree with Bootstrap Values")

# Add branch lengths (if needed)

p <- p + geom_treescale()

# Add bootstrap values to the branches (if they exist)
if (!is.null(rooted_tree$node.label) && length(rooted_tree$node.label) > 0 && !any(is.na(as.numeric(gsub("/.*", "", rooted_tree$node.label))))) {
  p <- p + geom_nodelab(aes(label = gsub("/.*", "", node.label)), size = 3, nudge_x = -0.01, nudge_y = 0.02, angle = 45)  # Adjust nudge and angle
}

p

dev.off()

print("Done!")