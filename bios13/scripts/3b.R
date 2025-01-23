# This script generates an offspring genome through recombination of two parent genomes.

# Clear the workspace.
rm(list = ls())

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

# Set a seed for the random number generator to ensure reproducible results.
set.seed(42)
parent1 = c(0, 1, 1, 0, 1) # Parent 1 genome
parent2 = c(1, 0, 0, 1, 0) # Parent 2 genome

offspring = recombination_func(parent1, parent2) # Call the above function here and assign the return of that function to offspring variable.
cat("Parent 1 Genome: ", parent1, "\n") # Display parent 1's genome.
cat("Parent 2 Genome: ", parent2, "\n") # Display parent 2's genome.
cat("Offspring Genome: ", offspring, "\n") # Display the resulting offspring genome.