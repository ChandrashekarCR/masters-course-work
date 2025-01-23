# This script simulates the differential equations metioned in the previous questions 4b and 4a, using a set of different initial conditions.
# i) The first plot is the concentration of X and Y with respect to time.
# ii) The second plot is in the phase plane together with the solutions in 4b.


# Clear the environment
rm(list = ls())

# Import Required Libraries
if(!require(ggplot2)){install.packages("ggplot2")}
if(!require(deSolve)){install.packages("deSolve")}
library(deSolve) # For solving the ODE
library(ggplot2)

# Define the ODE system for reaction dynamics
reaction_dynamics = function(t, state, parameters) {
  x = state[1] # Initial concentration of X (reactant)
  y = state[2] # Initial concentration of Y (product)
  
  k = parameters$k  # Reaction rate constant (positive)
  mu = parameters$mu  # Reaction rate constant (positive)
  
  dxdt = -(2 * k * x^2) + (2 * mu * y) # Rate of change of concentration of X (ODE 1)
  dydt = (k * x^2) - (mu * y) # Rate of change of concentration of Y (ODE 2)
  
  return(list(c(dxdt, dydt))) # deSolve package requires this particular format
}

# Define parameters and initial conditions
params = list(k = 0.5, mu = 0.2) # Any value can be given here
initial_conditions = list(
  c(x = 0.5, y = 0.1),
  c(x = 2, y = 1),
  c(x = 1, y = 2)  
)
time = seq(0, 10, by = 0.1) # Sequence of time points from 0 to 10 with increments of 0.1

# Solve the ODEs for each initial condition
solve_reaction_dynamics = function(initial_cond, params, time) {
  as.data.frame(ode(
    y = initial_cond,
    times = time,
    func = reaction_dynamics, # Pass the function defined previously here.
    parms = params
  ))
}

# Simulate dynamics for each initial condition and add labels
simulation_results = lapply(seq_along(initial_conditions), function(i) {
  # Solve the ODEs for the i-th initial condition
  result = solve_reaction_dynamics(initial_conditions[[i]], params, time)
  # Add a label indicating the initial condition for clarity
  result$Condition = paste("Condition", i)
  return(result)
})

# Combine simulation results from each of the initial conditions into a single data frame for easier plotting using ggplot2.
results_df = do.call(rbind, simulation_results)

# Generate caption text listing the initial conditions used
initial_conditions_caption = paste(
  "Initial Conditions:",
  paste(
    lapply(seq_along(initial_conditions), function(i) {
      # Format the initial conditions for each condition (X and Y values)
      sprintf("Condition %d: (x = %.1f, y = %.1f)", i, initial_conditions[[i]]["x"], initial_conditions[[i]]["y"])
    }),
    collapse = "; "  # Separate conditions with a semicolon and space
  )
)

# 4c i) Plot dynamics of X and Y over time
time_plot = ggplot(results_df, aes(x = time)) +
  geom_line(aes(y = x, color = "Concentration of X"), linewidth = 1) +
  geom_line(aes(y = y, color = "Concentration of Y"), linewidth = 1) +
  facet_wrap(~Condition, scales = "free_y") + # Doing a facet_wrap allows to plot all the three defined conditions side-by-side.
  labs(
    title = "Dynamics of X and Y Over Time",
    x = "Time",
    y = "Concentration",
    color = "Legend",
    caption = initial_conditions_caption # To know the initial values of each of the conditions, from the plot.
  ) +
  scale_color_manual(
    values = c(
      "Concentration of X" = "red",
      "Concentration of Y" = "blue"
    )
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),  # Centered and scaled title
    axis.title = element_text(size = 16, face = "bold"),               # Larger and bold axis titles
    axis.text = element_text(size = 14),                               # Larger axis text
    panel.grid.major = element_blank(),                                # Adjust grid line thickness
    panel.grid.minor = element_blank(),                                # Hide minor grid lines
    legend.text = element_text(size=12),                      
    legend.position = "bottom",  
    plot.caption = element_text(hjust = 0, face = "italic", size = 12)  # Caption styling
  )

# 4c ii) Plot 2: Phase plane with trajectories and equilibrium curve
plot_phase_plane <- function(simulations, params) {
  # Generate points for the equilibrium curve
  x_vals <- seq(0, 10, by = 0.1)
  y_vals <- (params$k / params$mu) * x_vals^2
  equilibrium_curve <- data.frame(x = x_vals, y = y_vals, label = "Equilibrium Line") # Add label for the legend
  
  # Combine simulation data from all conditions
  combined_data <- do.call(rbind, lapply(seq_along(simulations), function(i) {
    df <- simulations[[i]]
    df$Simulation <- paste("Condition", i)
    return(df)
  }))
  
  # Create the phase plane plot
  ggplot() +
    # Add the equilibrium curve with explicit color and label
    geom_line(data = equilibrium_curve, aes(x = x, y = y, color = label), 
              linetype = "dashed", linewidth = 1) +
    # Add trajectory paths for each initial concentration condition
    geom_path(data = combined_data, aes(x = x, y = y, color = Simulation), 
              linewidth = 1) +
    # Mark initial points
    geom_point(data = combined_data[combined_data$time == 0, ], 
               aes(x = x, y = y, color = Simulation), size = 3) +
    # Customize color scale to include equilibrium line and conditions
    scale_color_manual(
      name = "Legend",
      values = c(
        "Equilibrium Line" = "blue",
        "Condition 1" = "red",
        "Condition 2" = "green",
        "Condition 3" = "purple"
      )
    ) +
    labs(
      title = "Phase Plane with Trajectories and Equilibrium Curve",
      x = "Concentration of X",
      y = "Concentration of Y"
    ) +
    coord_cartesian(xlim = c(0, 2.5), ylim = c(0, 5)) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
      axis.title = element_text(size = 16, face = "bold"),
      axis.text = element_text(size = 14),
      legend.position = "bottom", # Position the legend at the bottom
      legend.title = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 12)
    )
}

# Generate the phase plane plot
phase_plane_plot = plot_phase_plane(simulation_results, params)

# Display plots
print(time_plot)
ggsave(plot = time_plot, filename = "./scripts/final_exam/exam_plots/time_plot.png", dpi = 600)

print(phase_plane_plot)
ggsave(plot = phase_plane_plot, filename = "./scripts/final_exam/exam_plots/phase_plane_traj.png", width = 200, units = "mm", dpi = 600)
