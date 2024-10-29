# Load necessary libraries
library(spatstat)
library(stelfi)
library(viridis)  
library(gridExtra)
library(sf)
library(reshape2)
library(INLA)
library(grid)
library(knitr)

# Parameters for Simulation
win <- owin(c(0, 5), c(0, 5)) 
beta0 <- 5
sigma2x <- 0.2  # Variance
range <- 1.2  # Spatial correlation range
nu <- 1
n_simulations <- 10  # Number of simulations

# Function to plot simulated point patterns
plot_simulated_point_patterns <- function(simulation, title) {
  lambda <- density(simulation, sigma = 0.1)  
  plot(lambda, main = title, col = inferno(128))  
  points(simulation, pch = 19, cex = 0.5, col = "black")  
}

# Generate PDF to plot the simulated point patterns for both INLA and STELFI models
pdf("Simulated_Point_Patterns_INLA_STELFI.pdf", width = 10, height = 6)
par(mfrow = c(1, 2))  

# Plot point pattern for INLA
lg.s_inla <- rLGCP('matern', beta0, var = sigma2x, scale = range / sqrt(8), nu = nu, win = win)
plot_simulated_point_patterns(lg.s_inla, "Simulated Point Pattern (INLA)")

# Plot point pattern for STELFI
lg.s_stelfi <- rLGCP('matern', beta0, var = sigma2x, scale = range / sqrt(8), nu = nu, win = win)
plot_simulated_point_patterns(lg.s_stelfi, "Simulated Point Pattern (STELFI)")

dev.off()
