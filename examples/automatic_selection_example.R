# Example: Using fit_movMFnoise for automatic model selection
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

# Use fit_movMFnoise to automatically select best G and noise option
fit_auto <- fit_movMFnoise(
  x,
  G = 1:4,
  noise = c(FALSE, TRUE),
  control = control_movMFnoise(nstart = 10, maxiter = 1000, noise_prop = 0.2),
  verbose = TRUE
)

cat("\n=== Best Model ===\n")
print(fit_auto)

cat("\n=== Summary ===\n")
print(summary(fit_auto))

# Compare with true labels
if (requireNamespace("mclust", quietly = TRUE)) {
  cat("\n=== Model Evaluation ===\n")
  ari <- mclust::adjustedRandIndex(
    true_labels,
    fit_auto$classification
  )
  cat(sprintf("Adjusted Rand Index: %.4f\n", ari))
}

cat("\n=== All Models Comparison ===\n")
print(fit_auto$model_summary)

cat("\n=== Visualize Model Selection ===\n")
summary_table <- fit_auto$model_summary

# Simple text-based visualization
cat("\nBIC by model:\n")
for (i in 1:nrow(summary_table)) {
  is_best <- (summary_table$G[i] == fit_auto$best_G &&
    summary_table$noise[i] == fit_auto$best_noise)

  noise_label <- ifelse(summary_table$noise[i], "+noise", "      ")
  marker <- ifelse(is_best, " <- BEST", "")
  converged_marker <- ifelse(summary_table$converged[i], "", " (not converged)")

  cat(sprintf(
    "  G=%d %s: BIC=%8.2f%s%s\n",
    summary_table$G[i],
    noise_label,
    summary_table$bic[i],
    marker,
    converged_marker
  ))
}

cat("\n=== Key Results ===\n")
cat(sprintf("Best number of components: %d\n", fit_auto$best_G))
cat(sprintf("Noise component included: %s\n", fit_auto$best_noise))
cat(sprintf("BIC: %.2f\n", fit_auto$bic))
cat(sprintf("ICL: %.2f\n", fit_auto$icl))
cat(sprintf("Log-likelihood: %.2f\n", fit_auto$loglik))

cat("\n=== Classification ===\n")
print(table(fit_auto$classification))
print(table(fit_auto$classification, true_labels))
cat("(0 = noise component)\n")
