#!/usr/bin/env Rscript
# Example: Fitting a Mixture of von Mises-Fisher Distributions with Noise

library(movMFnoise)

# Set random seed for reproducibility
set.seed(123)

cat("========================================\n")
cat("movMFnoise Example: Clustering on Sphere\n")
cat("========================================\n\n")

# Define true parameters for simulation
cat("Step 1: Define true model parameters\n")
cat("-------------------------------------\n")
theta_true <- list(
  alpha = c(0.3, 0.4),      # Mixing proportions
  mu = cbind(
    c(1, 0, 0),              # Mean direction 1
    c(0, 1, 0)               # Mean direction 2
  ),
  kappa = c(5, 10)          # Concentration parameters
)
cat("Components: 2\n")
cat("Noise proportion:", 1 - sum(theta_true$alpha), "\n\n")

# Generate synthetic data
cat("Step 2: Generate synthetic data\n")
cat("--------------------------------\n")
n <- 300
data <- rmovMFnoise(n = n, theta = theta_true)
cat("Generated", n, "observations in 3D\n")
cat("True cluster distribution:\n")
print(table(data$cluster))
cat("\n")

# Fit the model
cat("Step 3: Fit mixture model with EM algorithm\n")
cat("--------------------------------------------\n")
fit <- movMFnoise(
  x = data$data,
  G = 2,
  noise_init = 0.2,
  max_iter = 100,
  tol = 1e-6,
  verbose = FALSE
)

cat("Converged:", fit$converged, "\n")
cat("Iterations:", fit$iterations, "\n")
cat("Final log-likelihood:", round(fit$loglik[fit$iterations], 2), "\n\n")

# Display results
cat("Step 4: Compare estimated vs true parameters\n")
cat("---------------------------------------------\n")
cat("Mixing proportions (alpha):\n")
cat("  Estimated:", round(fit$theta$alpha, 3), "\n")
cat("  True:     ", theta_true$alpha, "\n\n")

cat("Noise proportion:\n")
cat("  Estimated:", round(fit$noise, 3), "\n")
cat("  True:     ", 1 - sum(theta_true$alpha), "\n\n")

cat("Concentration parameters (kappa):\n")
cat("  Estimated:", round(fit$theta$kappa, 2), "\n")
cat("  True:     ", theta_true$kappa, "\n\n")

cat("Estimated cluster distribution:\n")
print(table(fit$cluster))
cat("\n")

# Compute density at specific points
cat("Step 5: Evaluate density at test points\n")
cat("----------------------------------------\n")
test_points <- rbind(
  c(1, 0, 0),
  c(0, 1, 0),
  c(0, 0, 1)
)
rownames(test_points) <- c("Point 1", "Point 2", "Point 3")
dens <- dmovMFnoise(test_points, theta_true)

cat("Densities:\n")
for (i in 1:nrow(test_points)) {
  cat(" ", rownames(test_points)[i], ":", round(dens[i], 4), "\n")
}

cat("\n========================================\n")
cat("Example completed successfully!\n")
cat("========================================\n")
