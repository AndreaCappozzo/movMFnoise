# movMFnoise

R package for fitting Mixtures of von Mises-Fisher Distributions with a noise component (uniform density on the sphere).

## Overview

The `movMFnoise` package provides tools for modeling directional data on the unit sphere using mixtures of von Mises-Fisher distributions with an additional uniform noise component. This is particularly useful for:

- Clustering directional data with outliers
- Modeling spherical data with background noise
- Robust density estimation on the sphere

## Installation

You can install the development version from GitHub:

```r
# install.packages("devtools")
devtools::install_github("AndreaCappozzo/movMFnoise")
```

## Features

- **`rmovMFnoise()`**: Generate random samples from a mixture of von Mises-Fisher distributions with noise
- **`dmovMFnoise()`**: Compute density of the mixture model
- **`movMFnoise()`**: Fit the mixture model using the EM algorithm

## Usage

### Simulate data

```r
library(movMFnoise)

# Define model parameters
theta <- list(
  alpha = c(0.3, 0.4),  # Mixing proportions for vMF components
  mu = cbind(c(1, 0, 0), c(0, 1, 0)),  # Mean directions
  kappa = c(5, 10)  # Concentration parameters
)

# Generate 300 observations (30% noise)
set.seed(123)
data <- rmovMFnoise(n = 300, theta = theta)

# View cluster distribution
table(data$cluster)
```

### Fit a mixture model

```r
# Fit model with 2 components
fit <- movMFnoise(data$data, G = 2, noise_init = 0.1, verbose = TRUE)

# View estimated parameters
print(fit$theta)
print(fit$noise)

# View cluster assignments
table(fit$cluster)
```

### Compute density

```r
# Compute density for specific points
x <- rbind(c(1, 0, 0), c(0, 1, 0), c(0, 0, 1))
dens <- dmovMFnoise(x, theta)
print(dens)
```

## Model Details

The model assumes that observations come from a mixture of G von Mises-Fisher distributions plus a uniform noise component:

f(x) = Σ(g=1 to G) α_g * f_vMF(x | μ_g, κ_g) + α_noise * f_uniform(x)

where:
- α_g are the mixing proportions (sum to 1 - α_noise)
- f_vMF is the von Mises-Fisher density
- μ_g are the mean directions on the unit sphere
- κ_g are the concentration parameters
- f_uniform is the uniform density on the sphere

## License

GPL-3

## References

- Banerjee, A., Dhillon, I. S., Ghosh, J., & Sra, S. (2005). Clustering on the unit hypersphere using von Mises-Fisher distributions. Journal of Machine Learning Research, 6(Sep), 1345-1382.
