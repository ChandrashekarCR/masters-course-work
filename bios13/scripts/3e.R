# Clear environment
rm(list = ls())

# Importing Libraries
if(!require(ggplot2)){install.packages("ggplot2")}
library(ggplot2)
library(dplyr)

set.seed(15)

# Function for recombination
recombination_func = function(parent1, parent2) {
  if (length(parent1) != length(parent2)) {
    stop("Parent genomes must be of the same length")
  }
  
  genome_length = length(parent1)
  crossover_point = sample(1:(genome_length - 1), 1)
  offspring = c(parent1[1:crossover_point], parent2[(crossover_point + 1):genome_length])
  return(offspring)
}

# Simulation parameters
n = 100  # Population size
L = 100  # Genome length
generations = 1000

# Initialize the population
population = matrix(0, nrow = n, ncol = L)
population[1, ] = rep(1, L)  # First individual has a completely foreign genome

# Function to simulate one generation
simulate_generation = function(current_population, n, L) {
  next_population = matrix(0, nrow = n, ncol = L)
  
  for (i in 1:n) {
    parents_indices = sample(1:n, 2, replace = TRUE)
    parents = current_population[parents_indices, ]
    offspring = recombination_func(parents[1, ], parents[2, ])
    next_population[i, ] = offspring
  }
  
  return(next_population)
}

# Track the proportion of new genes and locus-specific frequencies
proportion_new_genes = numeric(generations)
locus_frequencies = matrix(0, nrow = generations, ncol = L)

# Run the simulation
for (gen in 1:generations) {
  proportion_new_genes[gen] = sum(population) / (n * L)
  locus_frequencies[gen, ] = colMeans(population)  # Track frequency at each locus
  population = simulate_generation(population, n, L)
}

# Data preparation for overall proportion plot
data = data.frame(
  Generation = 1:generations,
  Proportion = proportion_new_genes
)

# Plot the overall proportion of new genes
prop_gene_plot = ggplot(data, aes(x = Generation, y = Proportion)) +
  geom_line(color = "blue", linewidth = 1) +
  labs(
    title = "Proportion of New Genes Over Generations",
    x = "Generation",
    y = "Proportion"
  ) +
  theme_classic() + 
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),  # Centered and scaled title
    axis.title = element_text(size = 16, face = "bold"),               # Larger and bold axis titles
    axis.text = element_text(size = 14),                               # Larger axis text
    aspect.ratio = 0.6,                                                # Adjust aspect ratio
    panel.grid.major = element_blank(),                                # Hide major grid lines
    panel.grid.minor = element_blank()                                 # Hide minor grid lines
  )

# Prepare data for locus-specific frequencies plot (final generation)
final_gen_df = data.frame(
  Locus = 1:L,
  Frequency = locus_frequencies[generations, ]
)

# Plot locus-specific frequencies in the final generation
final_freq_plot = ggplot(final_gen_df, aes(x = Locus, y = Frequency)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  labs(
    title = "Frequency of New Genes at Each Locus (Final Generation)",
    x = "Locus",
    y = "Frequency"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),  # Centered and scaled title
    axis.title = element_text(size = 16, face = "bold"),               # Larger and bold axis titles
    axis.text = element_text(size = 14),                               # Larger axis text
    panel.grid.major = element_blank(),                                # Hide major grid lines
    panel.grid.minor = element_blank()                                 # Hide minor grid lines
  )

# Save the locus-specific plot
ggsave(plot = final_freq_plot, filename = "./scripts/final_exam/exam_plots/final_locus_frequencies.png", width = 200, dpi = 600, units = "mm")
ggsave(plot = prop_gene_plot, filename = "./scripts/final_exam/exam_plots/prop_gene_plot.png", width = 200, dpi = 600, units = "mm")

# Display the plots in RStudio
print(prop_gene_plot)
print(final_freq_plot)
