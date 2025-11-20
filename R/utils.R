#' Convert canonical parameters to mean direction and concentration
#'
#' @param theta Canonical parameter vector (theta = kappa * mu)
#' @return List containing:
#'   \item{mu}{Mean direction (unit vector)}
#'   \item{kappa}{Concentration parameter}
#' @export
canonical_to_mean <- function(theta) {
  kappa <- sqrt(sum(theta^2))
  mu <- if (kappa == 0) theta else theta / kappa
  return(list(mu = mu, kappa = kappa))
}

#' Compute MLE update for mu (mean direction)
#'
#' @param x Data matrix (n x d) where each row is a unit vector on the hypersphere
#' @return Normalized mean direction (unit vector)
#' @export
mu_update <- function(x) {
  # Compute the sample mean direction
  mean_direction <- colSums(x)

  # Normalize to unit length
  norm_val <- norm(mean_direction, type = "2")

  # Handle edge case where norm is zero
  mu <- if (norm_val > 0) {
    mean_direction / norm_val
  } else {
    mean_direction
  }

  return(mu)
}

#' Compute MLE update for kappa (concentration parameter)
#'
#' @param x Data matrix (n x d) where each row is a unit vector on the hypersphere
#' @param method Method for computing kappa (default: "Newton_Fourier")
#' @return Concentration parameter kappa
#' @export
kappa_update <- function(x, method = "Newton_Fourier") {
  n <- nrow(x)
  d <- ncol(x)

  # Compute the norm of the mean direction
  mean_direction <- colSums(x)
  R_bar <- norm(mean_direction, type = "2") / n

  # Get the solve_kappa function from movMF
  solve_kappa <- movMF:::get_solve_kappa(method)

  # Compute kappa using the specified method
  kappa <- solve_kappa$do_kappa(
    norms = R_bar * n,
    w = n,
    d = d,
    nu = 0
  )

  return(kappa)
}

#' Control Parameters for Mixture of von Mises-Fisher with Noise
#'
#' @description
#' Auxiliary function for controlling the EM algorithm used in fitting finite
#' mixtures of von Mises-Fisher distributions with an optional noise component.
#' This function specifies convergence criteria, initialization strategy, and
#' the method for computing concentration parameters.
#'
#' @param maxiter Integer. Maximum number of EM iterations. Default is 100.
#' @param reltol Numeric. Relative tolerance for convergence. The algorithm
#'   stops when the relative change in log-likelihood is smaller than \code{reltol}.
#'   Default is \code{sqrt(.Machine$double.eps)}.
#' @param kappa_method Character. Method for computing the concentration parameter
#'   kappa. Default is \code{"Newton_Fourier"}. See \code{\link[movMF]{movMF}}
#'   for available methods.
#' @param nstart Integer. Number of random starts for the EM algorithm. The best
#'   solution (highest log-likelihood) across all starts is returned. Default is 10.
#' @param start Character or matrix. Initialization method for the EM algorithm:
#'   \describe{
#'     \item{\code{"p"}}{Partition-based initialization (default)}
#'     \item{\code{"s"} or \code{"S"}}{Seeds-based initialization}
#'     \item{\code{"i"}}{Random class IDs initialization}
#'     \item{matrix}{User-provided matrix of initial posterior probabilities}
#'     \item{vector}{User-provided vector of initial class IDs}
#'   }
#'
#' @return A list with components corresponding to the control parameters.
#'
#' @seealso \code{\link{em_movMF}}
#'
#' @export
#'
#' @examples
#' # Default control parameters
#' control_movMFnoise()
#'
#' # Custom parameters for faster testing
#' control_movMFnoise(maxiter = 50, nstart = 5)
#'
#' # Tighter convergence tolerance with more iterations
#' control_movMFnoise(maxiter = 200, reltol = 1e-10)
#'
#' # Using seeds-based initialization
#' control_movMFnoise(start = "s")
control_movMFnoise <- function(
  maxiter = 100L,
  reltol = sqrt(.Machine$double.eps),
  kappa_method = "Newton_Fourier",
  nstart = 10L,
  start = "p",
  noise_prop = 0.1
) {
  list(
    maxiter = maxiter,
    reltol = reltol,
    kappa_method = kappa_method,
    nstart = nstart,
    start = start,
    noise_prop = noise_prop
  )
}
