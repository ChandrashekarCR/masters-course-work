# This script simulates the coupled dynamics of n1 and n2 and plots them as function of time.

# Clear the environment 
rm(list = ls())

# Import necessary libraries
# Check if libraries are installed, install if not
if (!require(ggplot2)) {
  install.packages("ggplot2")
}
if (!require(deSolve)) {
  install.packages("deSolve")
}
if (!require(tidyr)) {
  install.packages("tidyr")
}
library(deSolve)  # Load package for solving differential equations
library(ggplot2) # Load package for creating plots
library(tidyr)   # Load package for data manipulation

# Define a function for the coupled population dynamics
coupled_dynamics = function(t, state, params) {
  # Unpack the parameters from the 'params' list
  r1 = params["r1"] # Intrinsic growth of population 1
  r2 = params["r2"] # Intrinsic growth of population 2
  k = params["k"] # Strength of density dependence
  m = params["m"] # One directional dispersal rate from population 1 to population 2 (proportion)
  
  # Access state variables from the 'state' list
  n1 = state["n1"]
  n2 = state["n2"]
  
  # Define the population dynamics equations
  # dn1dt: rate of change of population 1
  # dn2dt: rate of change of population 2
  dn1dt = (r1 - k * n1^2) * n1 - m * n1
  dn2dt = (r2 - k * n2^2) * n2 + m * n1
  
  # Return a list containing the rates of change
  return(list(c(dn1dt, dn2dt)))
}

# Set the parameters for the model
params = c(r1 = 1, r2 = 0.6, k = 0.02, m = 0.1)

# Define the initial population sizes for each population
initial_state = c(n1 = 10, n2 = 3)

# Define the time steps for the simulation
timevec = seq(0, 100, by = 0.1)  # Sequence from 0 to 100 with steps of 0.1

# Solve the differential equations using the 'ode' function from deSolve
output = ode(y = initial_state, times = timevec, func = coupled_dynamics, parms = params)

# Convert the output from the 'ode' function to a data frame
output_df = as.data.frame(output)

# Reshape the data frame for easier plotting with ggplot2
# tidyr's pivot_longer creates separate columns for n1 and n2
population_data = pivot_longer(output_df, cols = c("n1", "n2"), names_to = "Population", values_to = "Size")

# Create the plot using ggplot2
coupled_dynamics_plot = ggplot(data = population_data, aes(x = time, y = Size, color = Population)) +
  scale_color_manual(values = c("n1" = "blue", "n2" = "red")) +  # Set colors for each population
  geom_line(linewidth = 1.5) +                                  # Plot lines for each population
  labs(
    title = "Coupled Dynamics of Two Populations (deSolve)",
    x = "Time",
    y = "Population Size",
    caption = paste(
      "Model: dn1/dt = (r1 - k*n1^2)n1 - m*n1",
      "\ndn2/dt = (r2 - k*n2^2)*n2 + m*n1",
      "\nParameters: r1 = 1, r2 = 0.6, k = 0.02, m = 0.1",
      "\nn1 =10, n2 = 3",
      sep = "\n"
    )  # Combine model equations, parameters, and initial conditions in caption
  ) +
  theme_classic() +                                               # Use a classic theme
  theme(
    plot.title = element_text(hjust = 0.5, vjust = 0.5, face = "bold", size = 18),  # Center and style title
    axis.title = element_text(face = "bold", size = 16),                         # Style axis titles
    plot.caption = element_text(hjust = 0, size = 12),                         # Style caption
    legend.text = element_text(size=12),                      # Hide minor grid lines
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size=14)
  )


# Save the plot as a high resolution image
ggsave(plot = coupled_dynamics_plot,filename = "./scripts/final_exam/exam_plots/coupled_dynamics.png", width = 200, units = "mm" , dpi=600)

# Display the plot
print(coupled_dynamics_plot)