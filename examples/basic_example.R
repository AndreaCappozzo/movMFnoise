# Example: Fitting a Mixture of von Mises-Fisher Distributions
# This example demonstrates how to use the em_movMF function

# Load required packages
library(movMFnoise)

# Set random seed for reproducibility
set.seed(123)

# Generate synthetic data on the unit sphere
# We'll create data from two distinct vMF components

n1 <- 50  # observations from component 1
n2 <- 50  # observations from component 2
d <- 3    # dimension (3D unit sphere)

# Component 1: concentrated around (1, 0, 0)
x1 <- matrix(rnorm(n1 * d), n1, d)
x1[, 1] <- x1[, 1] + 5  # Strong bias towards first dimension
x1 <- x1 / sqrt(rowSums(x1^2))  # Normalize to unit sphere

# Component 2: concentrated around (0, 1, 0)
x2 <- matrix(rnorm(n2 * d), n2, d)
x2[, 2] <- x2[, 2] + 5  # Strong bias towards second dimension
x2 <- x2 / sqrt(rowSums(x2^2))  # Normalize to unit sphere

# Combine data
x <- rbind(x1, x2)

# Fit mixture model with 2 components
cat("Fitting mixture of 2 von Mises-Fisher distributions...\n")
result <- em_movMF(x, G = 2, control = list(
  maxiter = 100,
  nstart = 5,
  verbose = TRUE
))

# Display results
cat("\n=== Results ===\n")
cat("Final log-likelihood:", result$loglik, "\n")
cat("Converged:", result$converged, "\n")
cat("Iterations:", result$iterations, "\n\n")

cat("Mixing proportions:\n")
print(result$parameters$pro)

cat("\nConcentration parameters (kappa):\n")
print(result$parameters$kappa)

cat("\nMean directions (mu):\n")
print(result$parameters$mu)

# Classify observations based on posterior probabilities
class_assignment <- apply(result$z, 1, which.max)
cat("\nClassification summary:\n")
print(table(class_assignment))

# Compare with true labels (first 50 are component 1, next 50 are component 2)
true_labels <- c(rep(1, n1), rep(2, n2))
cat("\nConfusion matrix:\n")
print(table(Predicted = class_assignment, True = true_labels))

cat("\nClassification accuracy:\n")
accuracy <- sum(class_assignment == true_labels | class_assignment == (3 - true_labels)) / length(true_labels)
print(accuracy)
