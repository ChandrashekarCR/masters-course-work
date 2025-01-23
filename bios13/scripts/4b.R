# This script plots the possible equilibrium in the xy phase plane.

# Clear the environment
rm(list = ls())

# Importing Libraries
if(!require(ggplot2)){install.packages("ggplot2")}
library(ggplot2)

# Define parameters
k = 0.5 # Rate constant
mu = 0.2 # Rate constant

# Generate a sequence of x values (non-negative as the chemical concentration cannot be less than zero.)
x_values = seq(0, 10, by = 0.1)

# Calculate corresponding y values using the equilibrium condition: y = (k/mu)*x^2
# Solve mathematically from previous question 4a.
y_values = (k / mu) * x_values^2

# Create a data frame for plotting
# ggplot2 requires the data to be in data frame format.
eqb_data = data.frame(x = x_values, y = y_values)

# Plot the equilibrium curve in the xy phase plane
eqb_curve = ggplot(data = eqb_data, aes(x, y)) + 
  geom_line(aes(color = "Equilibrium States"), linewidth = 1) +
  scale_color_manual(
    name = "Legend",
    values = c(
      "Equilibrium States" = "blue"
    )
  ) +
  labs(
    title = "Equilibrium States in the XY Phase Plane",
    x = "Concentration of X (x)",
    y = "Concentration of Y (y)"
  ) + 
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),  # Centered and scaled title
    axis.title = element_text(size = 16, face = "bold"),               # Larger and bold axis titles
    axis.text = element_text(size = 14),                               # Larger axis text
    panel.grid.major = element_blank(),                                 # Adjust grid line thickness
    panel.grid.minor = element_blank()                                  # Hide minor grid lines
  )

# Display the plot
print(eqb_curve)

# Save the plot as a high-resolution PNG
ggsave(
  filename = "./scripts/final_exam/exam_plots/eqb_phase_plane.png", 
  plot = eqb_curve, 
  width = 200,  # Width in mm for LaTeX integration
  units = "mm",
  dpi = 600
)
