#' Fit a Mixture of von Mises-Fisher Distributions with Noise
#'
#' @description
#' Fit a mixture of von Mises-Fisher distributions with a uniform noise component 
#' on the unit sphere using the EM algorithm.
#'
#' @param x Matrix of observations, each row is a point on the unit sphere.
#' @param G Integer. Number of von Mises-Fisher components.
#' @param noise_init Numeric. Initial proportion for the noise component (between 0 and 1).
#'                   Default is 0.1.
#' @param max_iter Integer. Maximum number of EM iterations. Default is 100.
#' @param tol Numeric. Convergence tolerance for the log-likelihood. Default is 1e-6.
#' @param init_method Character. Initialization method: "kmeans" or "random". 
#'                    Default is "kmeans".
#' @param verbose Logical. If TRUE, print iteration information. Default is FALSE.
#'
#' @return A list with components:
#' \itemize{
#'   \item{\code{theta}}{List of estimated parameters (alpha, mu, kappa)}
#'   \item{\code{noise}}{Estimated noise proportion}
#'   \item{\code{posterior}}{Matrix of posterior probabilities (n x (G+1))}
#'   \item{\code{cluster}}{Vector of hard cluster assignments (0 for noise, 1:G for components)}
#'   \item{\code{loglik}}{Vector of log-likelihood values at each iteration}
#'   \item{\code{iterations}}{Number of iterations performed}
#'   \item{\code{converged}}{Logical indicating whether the algorithm converged}
#' }
#'
#' @export
#'
#' @examples
#' # Generate data
#' theta_true <- list(
#'   alpha = c(0.3, 0.4),
#'   mu = cbind(c(1, 0, 0), c(0, 1, 0)),
#'   kappa = c(5, 10)
#' )
#' set.seed(123)
#' sim_data <- rmovMFnoise(n = 300, theta = theta_true)
#' 
#' # Fit model
#' fit <- movMFnoise(sim_data$data, G = 2)
#' print(fit$theta)
movMFnoise <- function(x, G, noise_init = 0.1, max_iter = 100, tol = 1e-6,
                       init_method = "kmeans", verbose = FALSE) {
  # Convert x to matrix if needed
  if (is.vector(x)) {
    x <- matrix(x, nrow = 1)
  }
  
  n <- nrow(x)
  d <- ncol(x)
  
  # Validate inputs
  if (G < 1) {
    stop("G must be at least 1")
  }
  if (noise_init < 0 || noise_init >= 1) {
    stop("noise_init must be between 0 and 1 (exclusive)")
  }
  
  # Initialize parameters
  init <- initialize_parameters(x, G, noise_init, init_method)
  alpha <- init$alpha
  mu <- init$mu
  kappa <- init$kappa
  noise <- init$noise
  
  # Initialize log-likelihood
  loglik <- rep(NA, max_iter)
  converged <- FALSE
  
  # EM algorithm
  for (iter in 1:max_iter) {
    # E-step: Compute posterior probabilities
    posterior <- e_step(x, alpha, mu, kappa, noise)
    
    # Compute log-likelihood
    loglik[iter] <- sum(log(rowSums(posterior)))
    
    if (verbose) {
      cat(sprintf("Iteration %d: loglik = %.4f\n", iter, loglik[iter]))
    }
    
    # Check convergence
    if (iter > 1) {
      if (abs(loglik[iter] - loglik[iter - 1]) < tol) {
        converged <- TRUE
        loglik <- loglik[1:iter]
        break
      }
    }
    
    # M-step: Update parameters
    m_result <- m_step(x, posterior)
    alpha <- m_result$alpha
    mu <- m_result$mu
    kappa <- m_result$kappa
    noise <- m_result$noise
  }
  
  # Hard cluster assignments
  cluster <- apply(posterior, 1, which.max) - 1  # 0 for noise, 1:G for components
  
  # Return results
  result <- list(
    theta = list(alpha = alpha, mu = mu, kappa = kappa),
    noise = noise,
    posterior = posterior,
    cluster = cluster,
    loglik = loglik,
    iterations = length(loglik),
    converged = converged
  )
  
  return(result)
}

#' Initialize Parameters for EM Algorithm
#' @keywords internal
initialize_parameters <- function(x, G, noise_init, init_method) {
  n <- nrow(x)
  d <- ncol(x)
  
  if (init_method == "kmeans") {
    # Use k-means for initialization
    km <- kmeans(x, centers = G, nstart = 10)
    cluster_assign <- km$cluster
  } else if (init_method == "random") {
    # Random initialization
    cluster_assign <- sample(1:G, size = n, replace = TRUE)
  } else {
    stop("init_method must be 'kmeans' or 'random'")
  }
  
  # Initialize mixing proportions
  alpha <- rep((1 - noise_init) / G, G)
  noise <- noise_init
  
  # Initialize mean directions
  mu <- matrix(0, nrow = d, ncol = G)
  for (g in 1:G) {
    cluster_data <- x[cluster_assign == g, , drop = FALSE]
    if (nrow(cluster_data) > 0) {
      mean_dir <- colMeans(cluster_data)
      mu[, g] <- mean_dir / sqrt(sum(mean_dir^2))
    } else {
      # Random direction if no points assigned
      mu[, g] <- rnorm(d)
      mu[, g] <- mu[, g] / sqrt(sum(mu[, g]^2))
    }
  }
  
  # Initialize concentration parameters
  kappa <- rep(1, G)
  for (g in 1:G) {
    cluster_data <- x[cluster_assign == g, , drop = FALSE]
    if (nrow(cluster_data) > 0) {
      r_bar <- sum(cluster_data %*% mu[, g]) / nrow(cluster_data)
      kappa[g] <- estimate_kappa(r_bar, d)
    }
  }
  
  return(list(alpha = alpha, mu = mu, kappa = kappa, noise = noise))
}

#' E-step: Compute Posterior Probabilities
#' @keywords internal
e_step <- function(x, alpha, mu, kappa, noise) {
  n <- nrow(x)
  G <- length(alpha)
  d <- ncol(x)
  
  # Compute component densities
  posterior <- matrix(0, nrow = n, ncol = G + 1)
  
  for (g in 1:G) {
    posterior[, g] <- alpha[g] * dmovMF_fast(x, mu[, g], kappa[g])
  }
  
  # Noise component (uniform on sphere)
  sphere_area <- 2 * pi^(d/2) / gamma(d/2)
  posterior[, G + 1] <- noise / sphere_area
  
  # Normalize
  posterior <- posterior / rowSums(posterior)
  
  return(posterior)
}

#' M-step: Update Parameters
#' @keywords internal
m_step <- function(x, posterior) {
  n <- nrow(x)
  G <- ncol(posterior) - 1
  d <- ncol(x)
  
  # Update mixing proportions
  n_g <- colSums(posterior)
  alpha <- n_g[1:G] / n
  noise <- n_g[G + 1] / n
  
  # Normalize alpha to ensure sum(alpha) + noise = 1
  total <- sum(alpha) + noise
  alpha <- alpha / total
  noise <- noise / total
  
  # Update mean directions and concentration parameters
  mu <- matrix(0, nrow = d, ncol = G)
  kappa <- rep(0, G)
  
  for (g in 1:G) {
    # Weighted mean direction
    weighted_sum <- colSums(posterior[, g] * x)
    r <- sqrt(sum(weighted_sum^2))
    
    if (r > 1e-10) {
      mu[, g] <- weighted_sum / r
      r_bar <- r / n_g[g]
      kappa[g] <- estimate_kappa(r_bar, d)
    } else {
      # If r is too small, keep previous values or use default
      mu[, g] <- rnorm(d)
      mu[, g] <- mu[, g] / sqrt(sum(mu[, g]^2))
      kappa[g] <- 1
    }
  }
  
  return(list(alpha = alpha, mu = mu, kappa = kappa, noise = noise))
}

#' Fast von Mises-Fisher Density Computation
#' @keywords internal
dmovMF_fast <- function(x, mu, kappa) {
  d <- length(mu)
  n <- nrow(x)
  
  # Normalizing constant
  C_d <- kappa^(d/2 - 1) / ((2*pi)^(d/2) * besselI(kappa, d/2 - 1))
  
  # Density
  inner_prod <- x %*% mu
  dens <- C_d * exp(kappa * inner_prod)
  
  return(as.vector(dens))
}

#' Estimate Kappa Using Approximation
#' @keywords internal
estimate_kappa <- function(r_bar, d) {
  # Approximation for kappa based on r_bar
  # Using the approximation: kappa ≈ (r_bar * d - r_bar^3) / (1 - r_bar^2)
  
  if (r_bar >= 1) {
    r_bar <- 0.999  # Avoid numerical issues
  }
  if (r_bar <= 0) {
    return(0.1)
  }
  
  # Simple approximation
  kappa <- r_bar * (d - r_bar^2) / (1 - r_bar^2)
  
  # Ensure kappa is positive and reasonable
  kappa <- max(0.1, min(kappa, 100))
  
  return(kappa)
}
