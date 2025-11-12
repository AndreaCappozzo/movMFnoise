#' Generate Random Samples from a Mixture of von Mises-Fisher Distributions with Noise
#'
#' @description
#' Generate random samples from a mixture of von Mises-Fisher distributions with 
#' a uniform noise component on the unit sphere.
#'
#' @param n Integer. Number of observations to generate.
#' @param theta List containing the mixture model parameters:
#' \itemize{
#'   \item{\code{alpha}}{Numeric vector of mixing proportions for the von Mises-Fisher 
#'         components (should sum to less than or equal to 1)}
#'   \item{\code{mu}}{Matrix where each column is a mean direction on the unit sphere 
#'         for each component}
#'   \item{\code{kappa}}{Numeric vector of concentration parameters (>= 0) for each component}
#' }
#' @param noise Numeric. Proportion of observations from the uniform noise component 
#'              (between 0 and 1). Default is 0, meaning \code{1 - sum(alpha)}.
#'
#' @return A list with components:
#' \itemize{
#'   \item{\code{data}}{Matrix of generated observations, each row is a point on 
#'         the unit sphere}
#'   \item{\code{cluster}}{Integer vector indicating the true cluster membership 
#'         (0 for noise, 1:G for components)}
#' }
#'
#' @export
#'
#' @examples
#' # Generate data from 2 components plus noise
#' theta <- list(
#'   alpha = c(0.3, 0.4),
#'   mu = cbind(c(1, 0, 0), c(0, 1, 0)),
#'   kappa = c(5, 10)
#' )
#' set.seed(123)
#' data <- rmovMFnoise(n = 300, theta = theta)
#' table(data$cluster)
rmovMFnoise <- function(n, theta, noise = NULL) {
  # Validate inputs
  if (!is.list(theta) || !all(c("alpha", "mu", "kappa") %in% names(theta))) {
    stop("theta must be a list with elements 'alpha', 'mu', and 'kappa'")
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
  
  # Validate dimensions
  if (length(kappa) != G) {
    stop("Length of kappa must match length of alpha")
  }
  if (ncol(mu) != G) {
    stop("Number of columns in mu must match length of alpha")
  }
  if (sum(alpha) > 1) {
    stop("Sum of alpha must be <= 1")
  }
  
  # Set noise proportion
  if (is.null(noise)) {
    noise <- 1 - sum(alpha)
  } else {
    if (noise < 0 || noise > 1) {
      stop("noise must be between 0 and 1")
    }
    if (abs(sum(alpha) + noise - 1) > 1e-10) {
      warning("sum(alpha) + noise should equal 1; adjusting alpha proportionally")
      alpha <- alpha * (1 - noise) / sum(alpha)
    }
  }
  
  # Generate cluster assignments
  probs <- c(alpha, noise)
  cluster <- sample(0:G, size = n, replace = TRUE, prob = probs)
  
  # Initialize data matrix
  data <- matrix(0, nrow = n, ncol = d)
  
  # Generate data from each component
  for (g in 1:G) {
    n_g <- sum(cluster == g)
    if (n_g > 0) {
      data[cluster == g, ] <- rmovMF(n_g, mu = mu[, g], kappa = kappa[g])
    }
  }
  
  # Generate noise observations (uniform on sphere)
  n_noise <- sum(cluster == 0)
  if (n_noise > 0) {
    # Generate uniform points on unit sphere
    noise_data <- matrix(rnorm(n_noise * d), nrow = n_noise, ncol = d)
    noise_data <- noise_data / sqrt(rowSums(noise_data^2))
    data[cluster == 0, ] <- noise_data
  }
  
  return(list(data = data, cluster = cluster))
}
