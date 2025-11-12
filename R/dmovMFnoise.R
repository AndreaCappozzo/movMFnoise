#' Density of a Mixture of von Mises-Fisher Distributions with Noise
#'
#' @description
#' Compute the density of a mixture of von Mises-Fisher distributions with 
#' a uniform noise component on the unit sphere.
#'
#' @param x Matrix of observations, each row is a point on the unit sphere.
#' @param theta List containing the mixture model parameters:
#' \itemize{
#'   \item{\code{alpha}}{Numeric vector of mixing proportions for the von Mises-Fisher 
#'         components}
#'   \item{\code{mu}}{Matrix where each column is a mean direction on the unit sphere 
#'         for each component}
#'   \item{\code{kappa}}{Numeric vector of concentration parameters for each component}
#' }
#' @param noise Numeric. Proportion of the uniform noise component. 
#'              Default is \code{1 - sum(alpha)}.
#' @param log Logical. If TRUE, return log-density. Default is FALSE.
#'
#' @return Numeric vector of density values for each observation.
#'
#' @export
#'
#' @examples
#' # Compute density for some points
#' theta <- list(
#'   alpha = c(0.4, 0.4),
#'   mu = cbind(c(1, 0, 0), c(0, 1, 0)),
#'   kappa = c(5, 10)
#' )
#' x <- rbind(c(1, 0, 0), c(0, 1, 0), c(0, 0, 1))
#' dmovMFnoise(x, theta)
dmovMFnoise <- function(x, theta, noise = NULL, log = FALSE) {
  # Validate inputs
  if (!is.list(theta) || !all(c("alpha", "mu", "kappa") %in% names(theta))) {
    stop("theta must be a list with elements 'alpha', 'mu', and 'kappa'")
  }
  
  # Convert x to matrix if needed
  if (is.vector(x)) {
    x <- matrix(x, nrow = 1)
  }
  
  alpha <- theta$alpha
  mu <- theta$mu
  kappa <- theta$kappa
  
  # Convert mu to matrix if it's a vector (single component case)
  if (is.vector(mu)) {
    mu <- matrix(mu, ncol = 1)
  }
  
  G <- length(alpha)  # Number of von Mises-Fisher components
  d <- nrow(mu)       # Dimension
  n <- nrow(x)        # Number of observations
  
  # Set noise proportion
  if (is.null(noise)) {
    noise <- 1 - sum(alpha)
  }
  
  # Compute density for each component
  dens <- matrix(0, nrow = n, ncol = G + 1)
  
  for (g in 1:G) {
    dens[, g] <- alpha[g] * dmovMF(x, mu = mu[, g], kappa = kappa[g])
  }
  
  # Add noise component (uniform density on d-dimensional sphere)
  # Surface area of unit sphere in d dimensions: 2*pi^(d/2) / gamma(d/2)
  sphere_area <- 2 * pi^(d/2) / gamma(d/2)
  dens[, G + 1] <- noise / sphere_area
  
  # Sum across components
  total_dens <- rowSums(dens)
  
  if (log) {
    return(log(total_dens))
  } else {
    return(total_dens)
  }
}
