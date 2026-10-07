# Example: Using movMFnoise for automatic model selection
# This demonstrates the high-level wrapper function for automatic model selection

# NOTE: If running this example from the package source directory, use:
#   devtools::load_all()
#   source("examples/automatic_selection_example.R")

library(movMFnoise)

# Set seed for reproducibility
set.seed(123)

# Simulate data on a 3D sphere
# 2 concentrated clusters + uniform noise

n_cluster1 <- 100
n_cluster2 <- 100
n_noise <- 30
n <- n_cluster1 + n_cluster2 + n_noise
d <- 3

# Cluster 1: concentrated around (1, 0, 0)
mu1 <- c(1, 0, 0)
kappa1 <- 10
data1 <- movMF::rmovMF(n_cluster1, theta = mu1 * kappa1)

# Cluster 2: concentrated around (-0.5, 0.866, 0) [120 degrees from cluster 1]
mu2 <- c(-0.5, 0.866, 0)
kappa2 <- 8
data2 <- movMF::rmovMF(n_cluster2, theta = mu2 * kappa2)

# Noise: uniform on sphere
data_noise <- matrix(rnorm(n_noise * d), ncol = d)
data_noise <- data_noise / sqrt(rowSums(data_noise^2))

# Combine all data
x <- rbind(data1, data2, data_noise)
true_labels <- c(rep(1, n_cluster1), rep(2, n_cluster2), rep(0, n_noise))

cat("=== Automatic Model Selection ===\n\n")

# Use movMFnoise to automatically select best G and noise option
fit_auto <- movMFnoise(data = x, G = 1:5,
                       noise = c(FALSE, TRUE),
                       control = control_movMFnoise(noise_prop = 0.2))

cat("\n=== Best Model ===\n")
print(fit_auto)

cat("\n=== Summary ===\n")
print(summary(fit_auto))

# Compare with true labels
if (requireNamespace("mclust", quietly = TRUE))
{
  cat("\n=== Model Evaluation ===\n")
  ari <- mclust::adjustedRandIndex(true_labels, fit_auto$classification)
  cat(sprintf("Adjusted Rand Index: %.4f\n", ari))
}

cat("\n=== Visualize Model Selection ===\n")

BIC = with(fit_auto$models_summary, as.data.frame(split(BIC, noise)))
matplot(BIC, type = "b", col = 1, pch = c(1, 19), lty = 1)
legend("bottomright", legend = c("Without Noise", "With Noise"),
       col = 1, pch = c(1, 19), lty = 1)

ICL = with(fit_auto$models_summary, as.data.frame(split(ICL, noise)))
matplot(ICL, type = "b", col = 1, pch = c(1, 19), lty = 1)
legend("bottomright", legend = c("Without Noise", "With Noise"),
       col = 1, pch = c(1, 19), lty = 1)
