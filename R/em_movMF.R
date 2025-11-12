# EM Algorithm for Mixture of vMF -----------------------------------------

#' EM Algorithm for Mixture of von Mises-Fisher Distributions
#'
#' @param x Data matrix (n x d) where each row is a unit vector on the hypersphere
#' @param G Number of mixture components (following mclust convention)
#' @param control List of control parameters:
#'   - maxiter: Maximum number of iterations (default: 100)
#'   - reltol: Relative tolerance for convergence (default: sqrt(.Machine$double.eps))
#'   - kappa_method: Method for computing kappa (default: "Banerjee_et_al_2005")
#'   - verbose: Print iteration details (default: FALSE)
#'   - nstart: Number of random starts (default: 10)
#' @return List with components:
#'   - parameters: List containing mu (G x d), kappa (length G), and pro (length G)
#'   - z: Posterior probabilities matrix (n x G)
#'   - loglik: Final log-likelihood value
#'   - control: List of control parameters used
#'   - n: Number of observations
#'   - d: Number of dimensions
#'   - G: Number of components
#' @export
em_movMF <- function(x, G, control = list()) {
  # Helper function for default values
  `%||%` <- function(x, y) if (is.null(x)) y else x

  # Normalize data to ensure unit vectors
  x <- x / sqrt(rowSums(x^2))

  n <- nrow(x)
  d <- ncol(x)

  # Set up control parameters with defaults
  control_default <- list(
    maxiter = 100L,
    reltol = sqrt(.Machine$double.eps),
    kappa_method = "Banerjee_et_al_2005",
    verbose = FALSE,
    nstart = 10L
  )

  # Override defaults with user-provided values
  control <- modifyList(control_default, control)

  # Extract control parameters
  maxiter <- control$maxiter
  reltol <- control$reltol
  kappa_method <- control$kappa_method
  verbose <- control$verbose
  nstart <- control$nstart

  # Get kappa solver
  solve_kappa_obj <- movMF:::get_solve_kappa(kappa_method)
  do_kappa <- solve_kappa_obj$do_kappa

  # Function to compute log-likelihood
  compute_loglik <- function(x, theta, pro) {
    # Compute cross_prod = x %*% t(theta) = <x_i, theta_g>
    cross_prod <- tcrossprod(x, theta)

    # Compute log normalizing constants
    # movMF uses: log f(x; theta) = <x, theta> - log(I_{d/2-1}(||theta||))
    # So log_C = -log(I_{d/2-1}(kappa))
    kappa_vals <- sqrt(rowSums(theta^2))
    log_C <- -log(besselI(kappa_vals, d / 2 - 1))

    # Log-likelihood contributions: log(pro_g) + <x_i, theta_g> + log_C_g
    log_densities <- sweep(cross_prod, 2, log(pro), "+")
    log_densities <- sweep(log_densities, 2, log_C, "+")

    # Sum over log-sum-exp for each observation
    max_vals <- apply(log_densities, 1, max)
    loglik <- sum(max_vals + log(rowSums(exp(log_densities - max_vals))))

    return(loglik)
  }

  # Function to initialize parameters
  initialize_params <- function() {
    # Random initialization: sample G observations as initial centers
    init_idx <- sample(n, G)
    mu_init <- x[init_idx, , drop = FALSE]
    mu_init <- mu_init / sqrt(rowSums(mu_init^2))

    # Initial kappa values (small random values)
    kappa_init <- runif(G, 1, 3)

    # Initial mixing proportions (uniform)
    pro_init <- rep(1 / G, G)

    # Construct theta = kappa * mu
    theta_init <- mu_init * kappa_init

    return(list(theta = theta_init, pro = pro_init))
  }

  # Run EM with multiple random starts
  best_result <- NULL
  best_loglik <- -Inf

  for (start in 1:nstart) {
    if (verbose) {
      cat(sprintf("\n=== Random start %d/%d ===\n", start, nstart))
    }

    # Initialize parameters
    init <- initialize_params()
    theta <- init$theta
    pro <- init$pro

    # EM iterations
    loglik_old <- -Inf
    converged <- FALSE

    for (iter in 1:maxiter) {
      # E-step: Compute posterior probabilities
      # cross_prod[i, g] = <x_i, theta_g>
      cross_prod <- tcrossprod(x, theta)

      # Compute log normalizing constants
      # movMF uses: log f(x; theta) = <x, theta> - log(I_{d/2-1}(||theta||))
      kappa_vals <- sqrt(rowSums(theta^2))
      log_C <- -log(besselI(kappa_vals, d / 2 - 1))

      # Compute log posteriors: log z[i,g] = log(pro_g) + <x_i, theta_g> + log_C_g
      log_z <- sweep(cross_prod, 2, log(pro), "+")
      log_z <- sweep(log_z, 2, log_C, "+")

      # Normalize to get probabilities (softmax)
      max_log_z <- apply(log_z, 1, max)
      z <- exp(log_z - max_log_z)
      z <- z / rowSums(z)

      # M-step: Update parameters

      # Update mixing proportions
      w <- colSums(z)
      pro <- w / n

      # Check for degenerate components (very small weight)
      if (any(pro < 1e-8)) {
        if (verbose) {
          cat("Warning: Degenerate component detected\n")
        }
        break
      }

      # Update mean directions
      # M[g, ] = sum_i z[i, g] * x[i, ]
      M <- crossprod(z, x) # G x d matrix

      # Compute norms
      norms <- sqrt(rowSums(M^2))

      # Normalize to get mu
      mu <- M / norms

      # Update kappa
      kappa <- do_kappa(norms, w, d, nu = 0)

      # Update theta = kappa * mu
      theta <- mu * kappa

      # Compute log-likelihood
      loglik <- compute_loglik(x, theta, pro)

      # Check convergence
      if (abs(loglik - loglik_old) / (abs(loglik) + reltol) < reltol) {
        converged <- TRUE
        if (verbose) {
          cat(sprintf("Converged at iteration %d\n", iter))
        }
        break
      }

      if (verbose && iter %% 10 == 0) {
        cat(sprintf("Iteration %d: loglik = %.4f\n", iter, loglik))
      }

      loglik_old <- loglik
    }

    # Store result if it's the best so far
    if (!is.na(loglik) && loglik > best_loglik) {
      best_loglik <- loglik
      best_result <- list(
        parameters = list(
          mu = mu,
          kappa = kappa,
          pro = pro
        ),
        z = z,
        loglik = loglik,
        control = control,
        n = n,
        d = d,
        G = G,
        iterations = iter,
        converged = converged
      )
    }
  }

  if (verbose) {
    cat(sprintf("\n=== Best result ===\n"))
    cat(sprintf("Final log-likelihood: %.4f\n", best_result$loglik))
  }

  class(best_result) <- "movMF"
  return(best_result)
}
