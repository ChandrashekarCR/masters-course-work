#!/usr/bin/env Rscript

rm(list=ls())

# List of necessary packages
packages <- c("ggplot2", "tidyverse", "dplyr")

# Install packages if they are not already installed
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

library(tidyverse)
library(dplyr)
library(ggplot2)

# Parse command-line arguments
args = commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript variants_based_statistics.R <input_prefix> <output_directory>")
}

input_prefix = args[1]  # Example: "./data/01_raw_data/subset_filter/subset_filter"
output_dir = args[2]  # Example: "./output"

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Function to save plot
save_plot <- function(plot, filename) {
  ggsave(filename, plot = plot, width = 8, height = 6, dpi = 300)
}

# Read and plot quality scores
var_qual = read_delim(paste0(input_prefix, ".lqual"), delim = "\t", col_names = c("chr","pos","qual"), skip=1)
var_qual_plot = ggplot(var_qual, aes(qual)) +
  geom_density(fill = "dodgerblue1", colour="black", alpha = 0.3) +
  theme_light()
save_plot(var_qual_plot, file.path(output_dir, "subset_filter.lqual.png"))

# Read and plot depth distribution
var_depth = read_delim(paste0(input_prefix, ".ldepth.mean"), delim = "\t", col_names = c("chr", "pos", "mean_depth", "var_depth"), skip = 1)
var_depth_plot = ggplot(var_depth, aes(mean_depth)) +
  geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) +
  theme_light()
save_plot(var_depth_plot, file.path(output_dir, "subset_filter.ldepth.mean.png"))

# Read and plot missingness
var_miss = read_delim(paste0(input_prefix, ".lmiss"), delim = "\t", col_names = c("chr", "pos", "nchr", "nfiltered", "nmiss", "fmiss"), skip = 1)
var_miss_plot = ggplot(var_miss, aes(fmiss)) + 
  geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) +
  theme_light()
save_plot(var_miss_plot, file.path(output_dir, "subset_filter.lmiss.png"))

# Read and plot individual depth
ind_depth = read_delim(paste0(input_prefix, ".idepth"), delim = "\t", col_names = c("ind", "nsites", "depth"), skip = 1)
ind_depth_plot = ggplot(ind_depth, aes(depth)) +
  geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3) +
  theme_light()
save_plot(ind_depth_plot, file.path(output_dir, "subset_filter.idepth.png"))

# Read and plot individual missingness
ind_miss = read_delim(paste0(input_prefix, ".imiss"), delim = "\t", col_names = c("ind", "ndata", "nfiltered", "nmiss", "fmiss"), skip = 1)
ind_miss_plot = ggplot(ind_miss, aes(fmiss)) +
  geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3) +
  theme_light()
save_plot(ind_miss_plot, file.path(output_dir, "subset_filter.imiss.png"))

# Read and plot heterozygosity
ind_het = read_delim(paste0(input_prefix, ".het"), delim = "\t", col_names = c("ind","ho", "he", "nsites", "f"), skip = 1)
ind_het_plot = ggplot(ind_het, aes(f)) + 
  geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3) + 
  theme_light()
save_plot(ind_het_plot, file.path(output_dir, "subset_filter.het.png"))

print("Done!")
