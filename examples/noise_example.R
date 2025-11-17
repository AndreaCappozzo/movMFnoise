# Example: Mixture of von Mises-Fisher with noise component
# This example demonstrates fitting a movMF model with a uniform noise component
# using fit_movMFnoise(), the main function of the package

# NOTE: If running this example from the package source directory, use:
#   devtools::load_all()
#   source("examples/noise_example.R")

library(movMFnoise)

# Set seed for reproducibility
set.seed(123)

# Simulate data on a 3D sphere
# 2 concentrated clusters + uniform noise

n_cluster1 <- 50
n_cluster2 <- 50
n_noise <- 30
n <- n_cluster1 + n_cluster2 + n_noise
d <- 3

# Cluster 1: concentrated around (1, 0, 0)
mu1 <- c(1, 0, 0)
kappa1 <- 10
theta1 <- kappa1 * mu1
data1 <- movMF::rmovMF(n_cluster1, theta = theta1)

# Cluster 2: concentrated around (-0.5, 0.866, 0) [120 degrees from cluster 1]
mu2 <- c(-0.5, 0.866, 0)
kappa2 <- 8
theta2 <- kappa2 * mu2
data2 <- movMF::rmovMF(n_cluster2, theta = theta2)

# Noise: uniform on sphere
# Generate random points and normalize
data_noise <- matrix(rnorm(n_noise * d), ncol = d)
data_noise <- data_noise / sqrt(rowSums(data_noise^2))

# Combine all data
x <- rbind(data1, data2, data_noise)
true_labels <- c(rep(1, n_cluster1), rep(2, n_cluster2), rep(0, n_noise))

# Fit model WITHOUT noise component
cat("=== Fitting model without noise ===\n")
fit_no_noise <- fit_movMFnoise(
  x,
  G = 2,
  noise = FALSE,
  control = control_movMFnoise(nstart = 5),
  verbose = TRUE,
)
print(fit_no_noise)

# Fit model WITH noise component
cat("\n=== Fitting model with noise ===\n")
fit_with_noise <- fit_movMFnoise(
  x,
  G = 2,
  noise = TRUE,
  control = control_movMFnoise(nstart = 5),
  verbose = TRUE
)
print(fit_with_noise)

# Compare classifications
cat("\n=== Classification comparison ===\n")
cat("True labels:\n")
print(table(true_labels))

cat("\nWithout noise component:\n")
print(table(fit_no_noise$classification))

cat("\nWith noise component:\n")
print(table(fit_with_noise$classification))

# Adjusted Rand Index (if package is available)
if (requireNamespace("mclust", quietly = TRUE)) {
  cat("\n=== Adjusted Rand Index ===\n")
  ari_no_noise <- mclust::adjustedRandIndex(
    true_labels,
    fit_no_noise$classification
  )
  ari_with_noise <- mclust::adjustedRandIndex(
    true_labels,
    fit_with_noise$classification
  )

  cat(sprintf("Without noise: %.4f\n", ari_no_noise))
  cat(sprintf("With noise: %.4f\n", ari_with_noise))
}

# Summary
cat("\n=== Model summaries ===\n")
cat("\nWithout noise:\n")
print(summary(fit_no_noise))

cat("\nWith noise:\n")
print(summary(fit_with_noise))
