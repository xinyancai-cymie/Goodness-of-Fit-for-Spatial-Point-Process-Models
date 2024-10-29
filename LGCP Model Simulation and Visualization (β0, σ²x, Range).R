# Load necessary libraries
library(spatstat)
library(ggplot2)
library(gridExtra)
library(INLA)

# Defining the Spatial Domain and Other Model Parameters
win <- owin(c(0, 5), c(0, 5))  # Define the spatial window

# Function to simulate LGCP model and generate the plot
simulate_and_plot <- function(beta0, sigma2x, range) {
  
  # Simulating an LGCP model with Matern covariance
  lg.s <- rLGCP('matern', beta0, var = sigma2x, scale = range / sqrt(8), nu = 1, win = win)
  
  # Create a pixel image from the intensity function
  intens <- density.ppp(lg.s, sigma = sqrt(sigma2x))
  
  # Convert the intensity to a dataframe
  intensity_data <- as.data.frame(intens)
  
  # Create a dataframe for the points
  points_data <- data.frame(x = lg.s$x, y = lg.s$y)
  
  # Generate the plot
  ggplot() +
    geom_raster(data = intensity_data, aes(x = x, y = y, fill = value)) +
    scale_fill_viridis_c(option = "C", direction = -1) +  # Using viridis color palette for better aesthetics
    geom_point(data = points_data, aes(x = x, y = y), color = "black", size = 2, shape = 1) +
    labs(title = paste("β0 =", beta0, "; σ²x =", sigma2x, "; range =", range),
         fill = "Intensity") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 15),
          axis.title.x = element_text(size = 12),
          axis.title.y = element_text(size = 12),
          panel.grid = element_blank())  # This removes the gray grid lines
}

# Set up the PDF output
pdf("LGCP_Simulated_Patterns_No_Grid.pdf", width = 10, height = 5)

# Plot for range = 0.3
p1 <- simulate_and_plot(beta0 = 6.0, sigma2x = 1.6, range = 0.3)

# Plot for range = 0.5
p2 <- simulate_and_plot(beta0 = 3.0, sigma2x = 0.02, range = 0.9)

# Arrange both plots in one figure
grid.arrange(p1, p2, ncol = 2)

# Close the PDF device
dev.off()


