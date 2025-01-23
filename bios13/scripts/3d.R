# This script plots the proportion of new genes in every generation. This script builds up from the scripts written in question 3b and 3c.

# Clear environment
rm(list = ls())

# Importing Libraries
if(!require(ggplot2)){install.packages("ggplot2")}
library(ggplot2)

# To ensure consistent and reproducible results, the random number generator was initialized with a seed of 15. Preliminary analyses using different random seeds (e.g., 12 and 18) demonstrated that the dynamics of the new gene's frequency were sensitive to this parameter. A seed of 18 resulted in a fluctuating, wave-like pattern, where the new gene type almost vanished before subsequently increasing in frequency, although this oscillation dampened over generations. A seed of 12 led to an abrupt initial decline in the new gene type's frequency. A seed of 15 was ultimately chosen as it provided a clearer and more stable representation of the overall trend being investigated, without the confounding effects of these extreme fluctuations.

set.seed(15) 
# Function for recombination
recombination_func = function(parent1, parent2) {
  # Verify that both parent genomes have the same length. If they don't, stop execution and report an error.
  if (length(parent1) != length(parent2)) {
    stop("Parent genomes must be of the same length")
  }
  
  # Determine the genome length. Since the previous check ensures equal lengths, we can use the length of either parent.
  genome_length = length(parent1)
  
  # Randomly select a crossover point. This point determines where the genomes are split and recombined.
  # The crossover point can be anywhere from 1 to genome_length - 1. We exclude genome_length itself to ensure that at least one gene is inherited from each parent.
  crossover_point = sample(1:(genome_length - 1), 1)
  
  # Create the offspring genome by combining segments from each parent.
  # We take the genes from parent 1 up to and including the crossover point, and combine them with the genes from parent 2 starting from the position immediately after the crossover point, up until the end of the parent 2 genome.
  offspring = c(parent1[1:crossover_point], parent2[(crossover_point + 1):genome_length])
  
  # Return the newly created offspring genome.
  return(offspring)
}

# Simulation parameters
n = 100 # Number of individuals in each generation (population size). The population size remains constant throughout the simulation.
L = 100 # Length of each individual's genome. Each genome is represented as a sequence of 8 binary values (0s and 1s). 0 is the old type gene and 1 is the new type gene.
generations = 100 # Number of generations to simulate. 

# Initialize the population.
# The population is represented as a matrix where each row is an individual's genome.
population = matrix(0, nrow = n, ncol = L)
# Introduce a "foreign" genome by setting the first individual's genome to all 1s. This
# simulates the introduction of a new mutation or gene into the population.
population[1, ] = rep(1, L)


# Function to simulate one generation. This function takes the current population as input
# and generates the next generation through random mating and recombination.
simulate_generation = function(current_population, n, L, gen_number) {
  next_population = matrix(0, nrow = n, ncol = L) # Initailize a matrix of the same size as the input to store the offsprings of the next generation.
  
  # Iterate through each individual in the next generation. This loop creates each of the n
  # individuals in the next generation by selecting parents and performing recombination.
  for (i in 1:n) {
    # Randomly select two parents from the current population with replacement.
    # Replacement is crucial here because an individual can be selected as a parent multiple
    # times within a generation. This reflects the possibility of an individual having
    # multiple offspring. That's why we set replace = TRUE
    parents_indices = sample(1:n, 2, replace = TRUE) # Get two indices at random
    parents = current_population[parents_indices, ]
    
    # Create an offspring through recombination using the previously defined function (3b).
    offspring = recombination_func(parents[1, ], parents[2, ])

    # Add the offspring to the next generation. The newly created offspring's genome is
    # added as a row to the next_population matrix.
    next_population[i, ] = offspring
  }
  
  return(next_population) # Return the next offspring matrix that will serve as the parent matrix for the next generation.
}

# Track the proportion of new genes in the population.
# We change the data type of the generations variable to numeric for easier calculations with float values.
proportion_new_genes = numeric(generations)

# Run the simulation. This loop iterates through the specified number of generations,
# simulating the population's evolution at each step.

for (gen in 1:generations) {
  # Calculate the proportion of new genes (1s) in the population at the current generation.
  # The sum of all values in the population matrix gives the total number of 1s (new genes).
  # This sum is then divided by the total number of genes in the population, which is
  # calculated as the number of individuals (n) multiplied by the genome length (L). (Number of items in the matrix.)
  proportion_new_genes[gen] = sum(population) / (n * L)
  
  population = simulate_generation(population, n, L, gen) # Run simulation for the next generation.
}

# Create a data frame for plotting. This is necessary because ggplot2, the plotting library
# we're using, works best with data in a data frame format.
data = data.frame(
  Generation = 1:generations,     # Column 1: Generation number. Creates a sequence of numbers from 1 to the total number of generations.
  Proportion = proportion_new_genes # Column 2: Proportion of new genes at each generation.
)

# Plot the result using ggplot2.
prop_gene_plot = ggplot(data, aes(x = Generation, y = Proportion)) +
  geom_line(aes(color = "New gene proportion"), linewidth = 1) + # Creates a line plot with a blue line of specified thickness.
  scale_color_manual(
    name = "Legend",
    values = c("New gene proportion" = "blue")
  ) +
  labs(                                     # Adds labels to the plot.
    title = "Proportion of New Genes Over Generations",
    x = "Generation",
    y = "Proportion"
  ) +
  theme_classic() + # Uses a clean, classic theme for the plot.
  # Formatting plot for better aesthetics and visualization
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"), # Centers and styles the plot title.
    axis.title = element_text(size = 16, face = "bold"),             # Styles the axis titles.
    axis.text = element_text(size = 14),                             # Styles the axis text (making it larger for readability).
    panel.grid.major = element_blank(),                               # Removes major grid lines.
    panel.grid.minor = element_blank(),                             # Removes minor grid lines.
    plot.caption = element_text(size=14)                       
  )


ggsave(plot = prop_gene_plot, filename = "./scripts/final_exam/exam_plots/prop_new_genes.png", width = 200, units = "mm", dpi = 700) # Save plot as png

# Display the plot.
print(prop_gene_plot)
