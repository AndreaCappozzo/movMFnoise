# movMFnoise

R package for fitting Mixture of von Mises-Fisher (vMF) distributions on the unit hypersphere using the EM algorithm.

## Installation

```r
# Install from GitHub (requires devtools)
devtools::install_github("AndreaCappozzo/movMFnoise")
```

## Usage

```r
library(movMFnoise)

# Create sample data on unit sphere
set.seed(123)
n <- 100
d <- 3
x <- matrix(rnorm(n * d), n, d)
x <- x / sqrt(rowSums(x^2))  # Normalize to unit sphere

# Fit mixture model with 2 components
result <- em_movMF(x, G = 2)

# Access results
result$parameters$mu     # Mean directions (G x d matrix)
result$parameters$kappa  # Concentration parameters (length G)
result$parameters$pro    # Mixing proportions (length G)
result$z                 # Posterior probabilities (n x G)
result$loglik            # Final log-likelihood

# Custom control parameters
result <- em_movMF(x, G = 2, control = list(
  maxiter = 100,
  nstart = 10,
  verbose = TRUE
))
```

## Features

- EM algorithm for mixture of von Mises-Fisher distributions
- Multiple random starts for better convergence
- Numerically stable log-space computations
- Configurable convergence criteria and iteration limits
- Compatible with movMF package conventions

## License

Apache License 2.0
