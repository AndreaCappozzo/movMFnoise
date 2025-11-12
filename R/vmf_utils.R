#' von Mises-Fisher Utility Functions
#'
#' @description
#' Internal utility functions for von Mises-Fisher distribution computations.
#'
#' @keywords internal
#' @name vmf_utils
NULL

#' Generate Random Samples from von Mises-Fisher Distribution
#'
#' @param n Integer. Number of observations to generate.
#' @param mu Numeric vector. Mean direction on the unit sphere.
#' @param kappa Numeric. Concentration parameter (>= 0).
#'
#' @return Matrix of generated observations, each row is a point on the unit sphere.
#'
#' @keywords internal
rmovMF <- function(n, mu, kappa) {
  d <- length(mu)
  
  if (kappa < 1e-10) {
    # If kappa is very small, generate uniform on sphere
    x <- matrix(rnorm(n * d), nrow = n, ncol = d)
    x <- x / sqrt(rowSums(x^2))
    return(x)
  }
  
  # Use the method from Wood (1994)
  # Generate samples in a coordinate system where mu = (0, 0, ..., 0, 1)
  # Then rotate to align with the actual mu
  
  samples <- matrix(0, nrow = n, ncol = d)
  
  for (i in 1:n) {
    # Generate W from the appropriate distribution
    w <- generate_w(kappa, d)
    
    # Generate uniform direction in d-1 dimensions
    v <- rnorm(d - 1)
    v <- v / sqrt(sum(v^2))
    
    # Combine to get sample in standard orientation
    x <- c(sqrt(1 - w^2) * v, w)
    
    # Rotate to align with mu
    samples[i, ] <- rotate_to_mu(x, mu)
  }
  
  return(samples)
}

#' Generate W Component for von Mises-Fisher Sampling
#' @keywords internal
generate_w <- function(kappa, d) {
  # Generate W using rejection sampling
  b <- (d - 1) / (2 * kappa + sqrt(4 * kappa^2 + (d - 1)^2))
  x0 <- (1 - b) / (1 + b)
  c <- kappa * x0 + (d - 1) * log(1 - x0^2)
  
  accept <- FALSE
  while (!accept) {
    z <- rbeta(1, (d - 1) / 2, (d - 1) / 2)
    w <- (1 - (1 + b) * z) / (1 - (1 - b) * z)
    u <- runif(1)
    
    if (kappa * w + (d - 1) * log(1 - x0 * w) - c >= log(u)) {
      accept <- TRUE
    }
  }
  
  return(w)
}

#' Rotate Vector to Align with Target Direction
#' @keywords internal
rotate_to_mu <- function(x, mu) {
  d <- length(mu)
  
  # If mu is already close to (0, 0, ..., 0, 1), no rotation needed
  if (abs(mu[d] - 1) < 1e-10) {
    return(x)
  }
  
  # Use Householder reflection to rotate
  # Create a vector that when normalized gives mu
  mu_normalized <- mu / sqrt(sum(mu^2))
  
  # Standard basis vector e_d = (0, 0, ..., 0, 1)
  e_d <- rep(0, d)
  e_d[d] <- 1
  
  # Householder vector
  v <- e_d - mu_normalized
  v_norm <- sqrt(sum(v^2))
  
  if (v_norm < 1e-10) {
    # mu is already e_d
    return(x)
  }
  
  v <- v / v_norm
  
  # Apply Householder transformation: (I - 2vv^T)x
  result <- x - 2 * v * sum(v * x)
  
  return(result)
}

#' Compute von Mises-Fisher Density
#'
#' @param x Matrix of observations, each row is a point on the unit sphere.
#' @param mu Numeric vector. Mean direction on the unit sphere.
#' @param kappa Numeric. Concentration parameter (>= 0).
#' @param log Logical. If TRUE, return log-density. Default is FALSE.
#'
#' @return Numeric vector of density values for each observation.
#'
#' @keywords internal
dmovMF <- function(x, mu, kappa, log = FALSE) {
  # Convert x to matrix if needed
  if (is.vector(x)) {
    x <- matrix(x, nrow = 1)
  }
  
  d <- length(mu)
  n <- nrow(x)
  
  if (kappa < 1e-10) {
    # Uniform density on sphere
    sphere_area <- 2 * pi^(d/2) / gamma(d/2)
    dens <- rep(1 / sphere_area, n)
    if (log) {
      return(log(dens))
    } else {
      return(dens)
    }
  }
  
  # Normalizing constant: C_d(kappa) = kappa^(d/2-1) / ((2*pi)^(d/2) * I_(d/2-1)(kappa))
  # where I_nu is the modified Bessel function of the first kind
  nu <- d / 2 - 1
  
  # Log normalizing constant
  log_C <- nu * log(kappa) - (d / 2) * log(2 * pi) - log(besselI(kappa, nu))
  
  # Inner product: x * mu
  inner_prod <- x %*% mu
  
  # Log density
  log_dens <- log_C + kappa * inner_prod
  
  if (log) {
    return(as.vector(log_dens))
  } else {
    return(as.vector(exp(log_dens)))
  }
}
