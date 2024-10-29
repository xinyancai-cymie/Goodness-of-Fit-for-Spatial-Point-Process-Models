# Goodness-of-Fit-for-Spatial-Point-Process-Models
This dissertation explores methods for assessing the goodness-of-fit (GOF) of spatial point process models, with a focus on Log-Gaussian Cox Processes (LGCPs). Using the Integrated Nested Laplace Approximation (INLA) for simulations, it evaluates several GOF measures—such as Pearson residuals, the Pair Correlation Function (PCF), the K-function, and global metrics like DIC and WAIC—under clustered and non-clustered spatial patterns. The study highlights both the strengths and limitations of these measures, particularly in capturing multi-scale clustering and localized spatial dependencies.

Data and code for this research are provided below:

**Data**

- 📊 **Earthquake_Events_Data.csv** - Dataset containing information on earthquake events (location, magnitude, date) used for LGCP modeling and spatial analysis.

**Code**

All other files contain R code scripts for different aspects of LGCP modeling and analysis:

- **INLA Parameter Analysis for LGCP Models (CPO, DIC, WAIC, Log-Likelihood).R** - This script performs parameter analysis for LGCP models using INLA, calculating metrics such as CPO, DIC, WAIC, and Log-Likelihood to assess model fit.

- **INLA Parameter Sensitivity Analysis for LGCP Models (PCF, K-Function, Residual Analysis).R** - R code focused on analyzing spatial patterns in LGCP models, using Pair Correlation Function (PCF), K-function, and residuals to evaluate clustering sensitivity.

- **LGCP Model Simulation and Visualization (β0, σ²x, Range).R** - Script for simulating LGCP models with varying parameters (β0, σ²x, and range) and visualizing the spatial point patterns to understand clustering behavior.

- **Simulated_Point_Patterns_Comparison_INLA_STELFI.R** - Compares point pattern simulations generated with INLA and STELFI, providing insights into differences in spatial clustering representation between the two methods.
