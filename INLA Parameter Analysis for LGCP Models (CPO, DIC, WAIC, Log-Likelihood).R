library(spatstat)
library(INLA)
library(ggplot2)
library(gridExtra)
library(sf)
library(reshape2)
library(igraph)
library(grid)

# Defining the Spatial Domain and Other Model Parameters
win <- owin(c(0, 5), c(0, 5))  # Define the spatial window
npix <- 300  # Number of pixels
beta0 <- 5  # Model intercept
sigma2x <- 0.07  # Variance
range <- 0.7  # Spatial correlation range
nu <- 1  # Smoothness parameter of the Matern covariance function
n_simulations <- 100  # Number of simulations

# Creating a 2D Mesh Representing the Spatial Domain
mesh <- inla.mesh.2d(loc.domain = cbind(c(0, 5), c(0, 5)), max.edge = c(0.5, 1), offset = c(0.1, 0.5), cutoff = 0.05)

# Defining the SPDE Model Using the Matern Covariance Function
spde <- inla.spde2.pcmatern(
  mesh = mesh,
  prior.range = c(1.2, 0.5),  # Prior for the range parameter
  prior.sigma = c(1, 0.01)  # Prior adjusted for larger sigma2x values
)

# Function to Analyze Model Fit for Different Sigma, Range, and Beta Values
analyze_model_fit <- function(sigma2x_vals, range_vals, beta_vals) {
  results <- data.frame()
  
  for (sigma2x in sigma2x_vals) {
    for (range in range_vals) {
      for (beta0 in beta_vals) {
        lg.s <- rLGCP('matern', beta0, var = sigma2x, scale = range / sqrt(8), nu = nu, win = win)
        xy <- cbind(lg.s$x, lg.s$y)[, 2:1]
        A <- inla.spde.make.A(mesh, xy)
        n <- nrow(xy)
        x <- rnorm(n, mean = 0, sd = sqrt(sigma2x))
        lambda <- exp(beta0 + x)  
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
        
        sim_result <- inla(y ~ 0 + intercept + f(s, model = spde),
                           family = 'poisson', 
                           data = inla.stack.data(stack), 
                           control.compute = list(cpo = TRUE, dic = TRUE, waic = TRUE, config = TRUE), 
                           control.predictor = list(A = inla.stack.A(stack)), 
                           E = inla.stack.data(stack)$e)
        
        fitted_values <- sim_result$summary.fitted.values$mean[1:n]
        residuals <- (y_sim - fitted_values) / sqrt(fitted_values)
        
        log_likelihood <- sum(log(sim_result$cpo$cpo), na.rm = TRUE)
        
        result <- data.frame(
          Sigma = sigma2x,
          Range = range,
          Beta = beta0,
          CPO = mean(sim_result$cpo$cpo, na.rm = TRUE),  # Handle any NA values
          DIC = sim_result$dic$dic,
          WAIC = sim_result$waic$waic,
          LogLikelihood = log_likelihood,
          ResidualMean = mean(residuals),
          ResidualSD = sd(residuals)
        )
        
        results <- rbind(results, result)
      }
    }
  }
  return(results)
}

sigma2x_vals <- c(0.02, 0.4, 0.8, 1.2, 1.6)
range_vals <- c(0.3, 0.5, 0.7, 0.9) 
beta_vals <- c(3.0, 4.0, 5.0, 6.0)  

# Analyze Model Fit for Different Sigma, Range, and Beta Values
fit_results <- analyze_model_fit(sigma2x_vals, range_vals, beta_vals)

print(fit_results)

p_cpo_avg <- ggplot(fit_results, aes(x = Sigma, y = CPO, color = as.factor(Beta))) +
  stat_summary(aes(group = as.factor(Beta)), geom = "line", fun = mean, size = 1.2) +
  stat_summary(aes(group = as.factor(Beta)), geom = "ribbon", fun.data = mean_cl_normal, fill = "grey80", alpha = 0.3, color = NA) +
  labs(title = "Average CPO for Different Sigma, Range, and Beta Values",
       x = "Sigma",
       y = "CPO",
       color = "Beta") +  
  theme_minimal() +
  scale_color_manual(values = c("darkred", "darkblue", "darkgreen", "darkorange")) +
  theme(legend.position = c(0.8, 0.8), 
        legend.background = element_rect(fill = "white", color = "black", size = 0.5))

p_dic_avg <- ggplot(fit_results, aes(x = Sigma, y = DIC, color = as.factor(Beta))) +
  stat_summary(aes(group = as.factor(Beta)), geom = "line", fun = mean, size = 1.2) +
  stat_summary(aes(group = as.factor(Beta)), geom = "ribbon", fun.data = mean_cl_normal, fill = "grey80", alpha = 0.3, color = NA) +
  labs(title = "Average DIC for Different Sigma, Range, and Beta Values",
       x = "Sigma",
       y = "DIC",
       color = "Beta") +  
  theme_minimal() +
  scale_color_manual(values = c("darkred", "darkblue", "darkgreen", "darkorange")) +
  theme(legend.position = c(0.2, 0.8),
        legend.background = element_rect(fill = "white", color = "black", size = 0.5)) 

p_waic_avg <- ggplot(fit_results, aes(x = Sigma, y = WAIC, color = as.factor(Beta))) +
  stat_summary(aes(group = as.factor(Beta)), geom = "line", fun = mean, size = 1.2) +
  stat_summary(aes(group = as.factor(Beta)), geom = "ribbon", fun.data = mean_cl_normal, fill = "grey80", alpha = 0.3, color = NA) +
  labs(title = "Average WAIC for Different Sigma, Range, and Beta Values",
       x = "Sigma",
       y = "WAIC",
       color = "Beta") +  
  theme_minimal() +
  scale_color_manual(values = c("darkred", "darkblue", "darkgreen", "darkorange")) +
  theme(legend.position = c(0.2, 0.8), 
        legend.background = element_rect(fill = "white", color = "black", size = 0.5)) 

p_loglik_avg <- ggplot(fit_results, aes(x = Sigma, y = LogLikelihood, color = as.factor(Beta))) +
  stat_summary(aes(group = as.factor(Beta)), geom = "line", fun = mean, size = 1.2) +
  stat_summary(aes(group = as.factor(Beta)), geom = "ribbon", fun.data = mean_cl_normal, fill = "grey80", alpha = 0.3, color = NA) +
  labs(title = "Average LogLikelihood for Different Sigma, Range, and Beta Values",
       x = "Sigma",
       y = "LogLikelihood",
       color = "Beta") +  # Changed legend label to "Beta"
  theme_minimal() +
  scale_color_manual(values = c("darkred", "darkblue", "darkgreen", "darkorange")) +
  theme(legend.position = c(0.9, 0.2), 
        legend.background = element_rect(fill = "white", color = "black", size = 0.5)) 

# Save the updated plots with Beta color distinctions and no external legend
pdf("INLA_Parameter_Analysis_Plots_Beta.pdf", width = 14, height = 10)

grid.arrange(
  grobs = list(
    p_cpo_avg,
    p_dic_avg,
    p_waic_avg,
    p_loglik_avg
  ),
  ncol = 2
)

dev.off()
