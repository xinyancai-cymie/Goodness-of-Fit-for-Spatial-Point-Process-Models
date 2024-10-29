# Load necessary libraries
library(spatstat)
library(INLA)
library(ggplot2)
library(gridExtra)
library(grid)

# Defining the Spatial Domain and Other Model Parameters
win <- owin(c(0, 5), c(0, 5))  # Define the spatial window
mesh <- inla.mesh.2d(loc.domain = cbind(c(0, 5), c(0, 5)), max.edge = c(0.5, 1), offset = c(0.1, 0.5), cutoff = 0.05)

# Defining the SPDE Model Using the Matern Covariance Function
spde <- inla.spde2.pcmatern(mesh = mesh, prior.range = c(1.2, 0.5), prior.sigma = c(1, 0.01))

# Selected parameter sets for spiky and non-spiky patterns
param_sets <- list(
  list(sigma2x = 1.6, range = 0.3, beta0 = 6.0),  # Spiky pattern
  list(sigma2x = 0.02, range = 0.9, beta0 = 3.0)  # Non-spiky pattern
)

# Function to analyze the selected parameter sets
analyze_specific_param_sets <- function(param_set, mesh, spde, win) {
  sigma2x <- param_set$sigma2x
  range <- param_set$range
  beta0 <- param_set$beta0
  
  # Simulating an LGCP model with Matern covariance
  lg.s <- rLGCP('matern', beta0, var = sigma2x, scale = range / sqrt(8), nu = 1, win = win)
  xy <- cbind(lg.s$x, lg.s$y)[, 2:1]
  
  # INLA model setup
  A <- inla.spde.make.A(mesh, xy)
  n <- nrow(xy)
  lambda <- exp(beta0 + rnorm(n, mean = 0, sd = sqrt(sigma2x)))
  y_sim <- rpois(n, lambda)
  
  y.pp <- c(rep(0, mesh$n), y_sim)
  e.pp <- rep(1, mesh$n + n)
  A.pp <- rbind(Diagonal(mesh$n, rep(1, mesh$n)), A)
  
  stack <- inla.stack(
    data = list(y = y.pp, e = e.pp),
    A = list(1, A.pp),
    effects = list(
      list(intercept = rep(1, mesh$n + n)),
      list(s = 1:spde$n.spde)
    ),
    tag = 'pp'
  )
  
  # Running the INLA model
  sim_result <- inla(y ~ 0 + intercept + f(s, model = spde),
                     family = 'poisson',
                     data = inla.stack.data(stack),
                     control.compute = list(cpo = TRUE, dic = TRUE, waic = TRUE, config = TRUE),
                     control.predictor = list(A = inla.stack.A(stack)),
                     E = inla.stack.data(stack)$e)
  
  # Return the result and simulated points for PCF and residuals
  return(list(result = sim_result, points = xy, lambda = lambda))
}

# Analyze spiky pattern
spiky_results <- analyze_specific_param_sets(param_sets[[1]], mesh, spde, win)

# Analyze non-spiky pattern
non_spiky_results <- analyze_specific_param_sets(param_sets[[2]], mesh, spde, win)

# Set up the PDF output 
pdf("Spiky_and_NonSpiky_Analysis.pdf", width = 10, height = 10)  

# Adjust margins and title spacing 
par(mfrow = c(1, 2), mar = c(6, 7, 5, 3), mgp = c(5, 2, 0), oma = c(2, 2, 3, 2))  

# Step 1: Pair Correlation Function (PCF) for spiky and non-spiky patterns
# Spiky Pattern
spiky_ppp <- ppp(spiky_results$points[, 1], spiky_results$points[, 2], window = win)
plot(pcf(spiky_ppp, rmax = 0.3), main = "Spiky Pattern - PCF", asp = 1)  # Setting aspect ratio to 1 for square plot

# Non-Spiky Pattern
non_spiky_ppp <- ppp(non_spiky_results$points[, 1], non_spiky_results$points[, 2], window = win)
plot(pcf(non_spiky_ppp, rmax = 0.3), main = "Non-Spiky Pattern - PCF", asp = 1)  # Setting aspect ratio to 1 for square plot

# Step 2: Residual analysis for spiky and non-spiky patterns
# Spiky residuals
residuals_spiky <- diagnose.ppm(ppm(spiky_ppp ~ 1), type = "pearson", plot.it = TRUE, legend = FALSE)
# Non-spiky residuals
residuals_non_spiky <- diagnose.ppm(ppm(non_spiky_ppp ~ 1), type = "pearson", plot.it = TRUE, legend = FALSE)

# Step 3: Envelope for K-function
# Spiky pattern
spiky_envelope <- envelope(spiky_ppp, fun = Kest, nsim = 99)
plot(spiky_envelope, main = "Spiky Pattern - K-function Envelope", asp = 1)

# Non-Spiky pattern
non_spiky_envelope <- envelope(non_spiky_ppp, fun = Kest, nsim = 99)
plot(non_spiky_envelope, main = "Non-Spiky Pattern - K-function Envelope", asp = 1)

# Close the PDF 
dev.off()
