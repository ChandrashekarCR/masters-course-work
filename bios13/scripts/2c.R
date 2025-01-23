# This script plots the rate of change of population as a function of population size n.

# Clear environment 
rm(list = ls())

# Load necessary libraries 
if(!require(ggplot2)){install.packages("ggplot2")}
library(ggplot2)

# Parameter Initialization
r = 1      # Intrinsic growth rate (rate of growth when resources are unlimited)
k = 0.02   # Strength of density dependence (how much the growth rate is affected by population density)

# Population size range (n)
# Creates a sequence of population sizes from 0 to 20 with increments of 0.1.
n = seq(0, 20, by = 0.1) 

# Population dynamics equation (defining a function is a good practice for reusability)
pop_eqn = function(n) {
  # This equation describes the rate of change of population (dn/dt) as a function of population size (n).
  return(r - (k * n^2)) * n
}

# Rate of change of population (dn/dt)
dndt = pop_eqn(n)

# Create a data frame for plotting (ggplot2 works best with data frames)
popdata = data.frame(n = n, dndt = dndt)
# Create a data frame for eqb line
eqb_line = data.frame(yintercept = 0)

# Create the plot using ggplot2
pop_growth_dynamics = ggplot(data = popdata, aes(x = n, y = dndt)) +
  geom_line(aes(color = "Population Dynamics"), linewidth = 1.2) + # Plots the line representing the population dynamics
  geom_hline(data = eqb_line, aes(yintercept = yintercept, color = "Equilibrium Points"), linetype = "dashed", linewidth = 1) + # Adds a horizontal red dashed line at y = 0. This line represents the equilibrium points where the population is not changing (dn/dt = 0).
  scale_color_manual(
    name = "Legend",
    values = c("Population Dynamics" = "blue",
               "Equilibrium Points" = "red")
  ) +
  labs(
    title = "Population Growth Dynamics",
    x = "Population Size (n)",
    y = "Rate of Change (dn/dt)",
    caption = paste("Model: dn/dt = (r - k * n^2) * n\n", # Using paste for better caption formatting
                    "Parameters: r =", r, ", k =", k) # Include parameter values in caption for clarity
  ) +
  theme_classic() + # Use a clean theme
  theme(
    plot.title = element_text(hjust = 0.5, vjust = 0.5, face = "bold", size = 18), # Center and style the title
    axis.title = element_text(face = "bold", size = 16), # Style axis titles
    axis.text = element_text(size = 14), # Style axis text
    legend.text = element_text(size = 12), # Style legend text
    plot.caption = element_text(hjust = 0, size = 12), # Style caption
    panel.grid.major = element_blank(), # Remove major grid lines for a cleaner look
    panel.grid.minor = element_blank()  # Remove minor grid lines
  )

# Save the plot (using relative path and specifying units and DPI is good practice)
ggsave(plot = pop_growth_dynamics, filename = "./scripts/final_exam/exam_plots/population_dynamics.png", width = 200, height = 150, units = "mm", dpi = 600) # Added height parameter to ggsave

# Print the plot to the console (useful for interactive sessions)
print(pop_growth_dynamics)