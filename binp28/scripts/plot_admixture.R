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

# Parse command-line arguments
args = commandArgs(trailingOnly = TRUE)
if (length(args) < 6) {
  stop("Usage: Rscript admixture_plot.R <file_prefix> <sample_path> <min_k> <max_k> <img_path> <individual_k>")
}

file_prefix = args[1]  # Example: "data/05_admixture/PQ_files/dataset_filter"
samplelist_path = args[2]
min_k = as.numeric(args[3])  # Example: 3
max_k = as.numeric(args[4])  # Example: 5
save_img = args[5]
individual_k = as.numeric(args[6])  # The specific k-value for individual plot

# Check for valid input
if (is.na(min_k) | is.na(max_k) | min_k > max_k) {
  stop("Invalid K range. Ensure min_k and max_k are numeric and min_k <= max_k.")
}

if (is.na(individual_k) | individual_k < min_k | individual_k > max_k) {
  stop("Invalid individual_k value. Ensure individual_k is between min_k and max_k.")
}

# Read sample list
if (!file.exists(samplelist_path)) {
  stop(paste("Sample list file not found:", samplelist_path))
}

samplelist = read_tsv(samplelist_path, col_names = "sample", show_col_types = FALSE)

# Initialize tibble to store all data
all_data = tibble(sample=character(), k=numeric(), Q=character(), value=numeric())

# Loop over K values to read data
for (k in seq(min_k, max_k)) {
  file_path = paste0(file_prefix, ".", k, ".Q")
  
  if (!file.exists(file_path)) {
    warning(paste("File not found, skipping:", file_path))
    next
  }
  
  # Read Q matrix
  data = read_delim(file_path, 
                    col_names = paste0("Q", seq(1:k)), 
                    delim = " ", 
                    show_col_types = FALSE)
  
  # Add sample and k values
  if (nrow(data) != nrow(samplelist)) {
    stop(paste("Mismatch between sample list and Q matrix for k =", k))
  }
  
  data$sample = samplelist$sample
  data$k = k
  
  # Reshape to long format
  data = data %>% pivot_longer(cols = starts_with("Q"), names_to = "Q", values_to = "value")
  
  # Append to all_data
  all_data = bind_rows(all_data, data)
}

# Plot all K values with facet_wrap
p_all = ggplot(all_data, aes(x = sample, y = value, fill = factor(Q))) + 
  geom_bar(stat = "identity", position = "stack") +
  xlab("Sample") + ylab("Ancestry") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 16, hjust = 0.5),
        legend.text = element_text(size = 14)
  )+
  scale_fill_brewer(palette = "Set1", name = "Cluster") +
  facet_wrap(~k, ncol = 1)

# Save the plot with all K values
ggsave(save_img, plot = p_all, width = 12, height = 6, dpi = 300, device = "png")


# Extract directory and base name from save_img
save_dir <- dirname(save_img)
save_base <- tools::file_path_sans_ext(basename(save_img))

# Plot for individual k value
p_individual = ggplot(all_data[all_data$k == individual_k, ], aes(x = sample, y = value, fill = factor(Q))) + 
  geom_bar(stat = "identity", position = "stack") +
  xlab("Sample") + ylab("Ancestry") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
          axis.text = element_text(size = 14),
          axis.title = element_text(size = 14),
          plot.title = element_text(size = 16, hjust = 0.5),
        legend.text = element_text(size = 14)
        )+
  scale_fill_brewer(palette = "Set1", name = "Cluster") +
  ggtitle(paste("Individual Plot for K =", individual_k))

# Construct the new path for the individual plot
individual_plot_path <- file.path(save_dir, paste0(save_base, "_individual_k_", individual_k, ".png"))

# Save the individual plot
ggsave(individual_plot_path, plot = p_individual, width = 12, height = 6, dpi = 300, device = "png")


print("All plots saved successfully!")



