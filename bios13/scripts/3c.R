# This script simulates the evolutionary dynamics of a population's genome over multiple generations,
# using a simplified model of genetic recombination. The simulation tracks the spread of a "foreign"
# gene within a population of fixed size.
# We can always change the generation to 100, but for clarity and conciseness of the material we will limit the generations to 10 
# and see the genomic population structure at every generation to see the old type and new type genes in the genome.


# Clear environment
rm(list = ls())

# Set the random state
set.seed(18)
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
n = 10 # Number of individuals in each generation (population size). The population size remains constant throughout the simulation.
L = 8  # Length of each individual's genome. Each genome is represented as a sequence of 8 binary values (0s and 1s). 0 is the old type gene and 1 is the new type gene.
generations = 10 # Number of generations to simulate. Reduced from a potential 100 for clearer output in this example. Running for 10 generations allows us to easily observe the changes in the population's genetic makeup.

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
  
  cat("Generation: ", gen_number, "\n")
  
  # Iterate through each individual in the next generation. This loop creates each of the n
  # individuals in the next generation by selecting parents and performing recombination.
  for (i in 1:n) {
    # Randomly select two parents from the current population with replacement.
    # Replacement is crucial here because an individual can be selected as a parent multiple
    # times within a generation. This reflects the possibility of an individual having
    # multiple offspring. That's why we set replace = TRUE
    parents_indices = sample(1:n, 2, replace = TRUE) # Get two indices at random
    parents = current_population[parents_indices, ]
    
    cat("Parents for offspring ", i, ":\n") # Display the parents that are considered for generating an offspring.
    cat("Parent 1 (Index ", parents_indices[1], "):", parents[1, ], "\n")
    cat("Parent 2 (Index ", parents_indices[2], "):", parents[2, ], "\n")
    
    # Create an offspring through recombination using the previously defined function (3b).
    offspring = recombination_func(parents[1, ], parents[2, ])
    cat("Offspring ", i, ":", offspring, "\n\n")
    
    # Add the offspring to the next generation. The newly created offspring's genome is
    # added as a row to the next_population matrix.
    next_population[i, ] = offspring
  }
  
  return(next_population) # Return the next offspring matrix that will serve as the parent matrix for the next generation.
}

# Run the simulation. This loop iterates through the specified number of generations,
# simulating the population's evolution at each step.
for (gen in 1:generations) {
  population = simulate_generation(population, n, L, gen)
  cat("Population at the end of Generation ", gen, ":\n") # Display the genomic population structure to see how the old type (0) and new type (1) genes are seen.
  print(population)
  cat("\n-------------------------------------------\n") # For better aesthetics and seperation for each generation.
}