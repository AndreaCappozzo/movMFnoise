#' Fit Mixture of von Mises-Fisher Distributions with Automatic Model Selection
#'
#' @description
#' Wrapper function that fits mixture of von Mises-Fisher distributions for
#' multiple values of G (number of components) and selects the best model
#' according to BIC, similar to \code{mclust::Mclust()}. The main functionality
#' includes the optional addition of a uniform noise component to handle
#' outliers and observations that do not belong to any cluster.
#'
#' @param data Data matrix (n x d) where each row is a unit vector on the
#'   hypersphere. If not, rows will be normalized to unit length.
#' @param G Vector of integers specifying the number of mixture components to try.
#'   If `NULL` (default) sets `G = 1:9`.
#' @param noise Logical or character vector; if `TRUE`, fits models with a
#'   uniform noise component. If `FALSE` (default), fits models without noise.
#'   If `c(TRUE, FALSE)`, tries both options for each G and selects the best
#'   overall model.
#' @param control List of control parameters. See \code{\link{control_movMFnoise}}
#'   for details. Default values are obtained by calling \code{control_movMFnoise()}.
#' @param verbose Logical; if `TRUE`, print progress information (default: `FALSE`)
#' @param criterion Character; model selection criterion, either `"BIC"` (default)
#'   or `"ICL"`.
#' @param keep_all Logical; if `TRUE` all fitted models are returned in the output
#'   (default: `FALSE`).
#'
#' @return An object of class `'movMFnoise'` containing the elements provided
#'   by \code{\link{em_movMF}} for the best model, plus additional information
#'   about all fitted models:
#'   \item{models_summary}{Data frame of summary statistics for model selection;}
#'   \item{models}{List of all fitted models (optional, only if `keep_all = TRUE`)}
#'
#' @details
#' This function fits movMF models for each combination of G values and noise
#' options, then selects the model with the best (highest) BIC or ICL value.
#' This follows the approach used by \code{mclust::Mclust()} for Gaussian mixtures,
#' where higher values of BIC and ICL indicate better models.
#'
#' The function returns the best model along with summary information about all
#' fitted models for comparison purposes.
#'
#' @seealso \code{\link{em_movMF}}, \code{\link{control_movMFnoise}}
#'
#' @examples
#' \dontrun{
#' # Generate example data
#' set.seed(123)
#' x <- movMF::rmovMF(100, theta = rbind(c(5, 0, 0), c(0, 5, 0)))
#'
#' # Automatic model selection
#' fit <- movMFnoise(x, G = 1:5)
#' print(fit)
#' summary(fit)
#'
#' # Try both with and without noise
#' fit <- movMFnoise(x, G = 1:5, noise = c(FALSE, TRUE))
#'
#' # Compare all models
#' plot(fit$G_sequence, fit$bic_all, type = "b", xlab = "G", ylab = "BIC")
#' }
#' @export

movMFnoise <- function(data, G = NULL, noise = FALSE,
                       control = control_movMFnoise(),
                       verbose = interactive(),
                       criterion = c("BIC", "ICL"),
                       keep_all = FALSE)
{
  x <- data.matrix(data)
  # normalize data to ensure unit vectors
  x <- x / sqrt(rowSums(x^2))

  # check G using default as in mclust
  G = if(is.null(G)) 1:9 else sort(as.integer(unique(G)))

  # check noise specification
  if (length(noise) == 1)
  {
    noise_options <- noise
  } else
  {
    # Try both options
    noise_options <- unique(noise)
  }

	# check criterion
  criterion <- match.arg(criterion, several.ok = FALSE)

  # Create grid of models to fit
  model_grid <- expand.grid(G = G,
                            noise = noise_options,
                            stringsAsFactors = FALSE)
  n_models <- nrow(model_grid)

  if(verbose)
    cat("Fitting movMF models... \n")

  # Initialize storage
  models <- vector("list", n_models)
  loglik_all <- numeric(n_models)
  df_all <- numeric(n_models)
  bic_all <- numeric(n_models)
  icl_all <- numeric(n_models)
  converged_all <- logical(n_models)

  # Fit all models
  for (i in 1:n_models)
  {
    G_i <- model_grid$G[i]
    noise_i <- model_grid$noise[i]

    if(verbose)
    {
      cat(sprintf("  Model %d/%d: G=%d, noise=%s ... ",
                  i, n_models, G_i, noise_i))
    }

    fit_i <- try(em_movMF(data = x,
                          G = G_i,
                          control = control,
                          noise = noise_i,
                          verbose = FALSE),
                 silent = TRUE)

    if(inherits(fit_i, "try-error"))
    {
      models[[i]] <- NA
      loglik_all[i] <- NA
      df_all[i] <- NA
      bic_all[i] <- NA
      icl_all[i] <- NA
      converged_all[i] <- FALSE
    } else
    {
      models[[i]] <- fit_i
      loglik_all[i] <- fit_i$loglik
      df_all[i] <- fit_i$df
      bic_all[i] <- fit_i$bic
      icl_all[i] <- fit_i$icl
      converged_all[i] <- fit_i$converged
      if(verbose)
      {
        cat(sprintf("BIC=%.*g, ICL=%.*g, converged=%s\n",
                    getOption("digits"), fit_i$bic,
                    getOption("digits"), fit_i$icl,
                    fit_i$converged))
      }
    }
  }

  # Select best model based on criterion
  criterion_values <- if(criterion == "BIC") bic_all else icl_all

  # Only consider converged models
  valid_models <- converged_all & !is.na(criterion_values)

  if(!any(valid_models))
    stop("No models converged successfully")

  # Find best model (highest BIC/ICL, following mclust convention)
  best_idx <- which.max(criterion_values[valid_models])
  # Map back to original index
  best_idx <- which(valid_models)[best_idx]

  best_model <- models[[best_idx]]
  best_G <- model_grid$G[best_idx]
  best_noise <- model_grid$noise[best_idx]

  if (verbose & n_models > 1)
  {
    cat(sprintf("Best model: G=%d, noise=%s, %s=%.2f\n",
                best_G, best_noise,
                toupper(criterion), criterion_values[best_idx]))
  }

  # Add model selection information to best model
  # best_model$G_sequence <- model_grid$G
  # best_model$noise_sequence <- model_grid$noise
  # best_model$bic_all <- bic_all
  # best_model$icl_all <- icl_all
  # best_model$loglik_all <- loglik_all
  # best_model$converged_all <- converged_all
  # best_model$best_G <- best_G
  # best_model$best_noise <- best_noise
  # best_model$noise <- best_noise
  best_model$criterion <- criterion

  # Add model selection summary table
  models_summary <- data.frame(G = model_grid$G,
                               noise = model_grid$noise,
                               loglik = loglik_all,
                               df = df_all,
                               BIC = bic_all,
                               ICL = icl_all,
                               converged = converged_all,
                               stringsAsFactors = FALSE)
  best_model$models_summary <- models_summary

  # Optionally include all fitted models
  if(keep_all)
    best_model$models <- models

  class(best_model) <- "movMFnoise"
  return(best_model)
}

#' EM Algorithm for Mixture of von Mises-Fisher Distributions with Noise
#'
#' @param data Data matrix (n x d) where each row is a unit vector on the hypersphere
#' @param G Number of mixture components (following mclust convention)
#' @param control List of control parameters. See \code{\link{control_movMFnoise}} for details.
#'   Default values are obtained by calling \code{control_movMFnoise()}.
#' @param noise Logical or numeric vector; if TRUE, a uniform noise component is included.
#'   If numeric, specifies which observations are initially assigned to noise (default: FALSE)
#' @param verbose Logical; if TRUE, print iteration details (default: FALSE)
#' @return List with components:
#'   \item{data}{Input data matrix.}
#'   \item{n}{Number of observations.}
#'   \item{d}{Number of dimensions.}
#'   \item{G}{Number of components (excluding noise).}
#'   \item{parameters}{List containing:
#'      `mu` (d x G matrix with named columns 1:G),
#'     `kappa` (named vector of length G with names 1:G),
#'     `pro` (named vector of length G+1 if noise, else G, with names 1:G or c(1:G, 0) if noise), and
#'     `Vinv` (inverse hypervolume of noise component, if `noise=TRUE`).}
#'   \item{hypvol}{Hypervolume (1/Vinv) if noise component is present, else `NA`.}
#'   \item{z}{Posterior probabilities matrix (n x (G+1) if noise, else n x G).}
#'   \item{classification}{Vector of MAP component assignments (1:G, or 0 for noise if present).}
#'   \item{loglik}{Final log-likelihood value.}
#'   \item{loglik_trace}{Vector of log-likelihood values at each EM iteration.}
#'   \item{bic}{Bayesian Information Criterion (higher is better, following mclust convention).}
#'   \item{icl}{Integrated Completed Likelihood criterion (higher is better).}
#'   \item{control}{List of control parameters used.}
#'   \item{iterations}{Number of iterations performed.}
#'   \item{converged}{Logical indicating convergence.}
#'   \item{criterion}{String specifying model selection criterion used.}
#'   \item{models_summary}{Data frame containing summary of fitted models.}
#' @importFrom utils modifyList
#' @importFrom stats runif
#' @export
em_movMF <- function(data, G, noise = FALSE,
                     control = control_movMFnoise(),
                     verbose = FALSE)
{
  x <- data.matrix(data)
  n <- nrow(x)
  d <- ncol(x)
  # normalize data to ensure unit vectors
  x <- x / sqrt(rowSums(x^2))

  # Handle noise specification
  has_noise <- FALSE
  noise_idx <- NULL
  if (is.logical(noise) && length(noise) == 1) {
    has_noise <- noise
  } else if (is.numeric(noise)) {
    has_noise <- TRUE
    noise_idx <- noise
  }

  # Compute inverse hypervolume (Vinv) for uniform noise on unit sphere
  # Surface area of unit (d-1)-sphere: S_{d-1} = 2*pi^(d/2) / Gamma(d/2)
  # Vinv = 1 / S_{d-1}
  # IMPORTANT: the normalized surface-area measure is already included in
  # the vMF specification of the movMF package, so we do not need to
  # multiply by any additional constants here and the uniform density is simply 1.
  Vinv <- if (has_noise) {
    # gamma(d / 2) / (2 * pi^(d / 2)) # FIXME double-check this constant
    1
  } else {
    NULL
  }

  # Merge user-provided control parameters with defaults
  # If control is already a complete list, use modifyList to override defaults
  control <- modifyList(control_movMFnoise(), control)

  # Extract control parameters
  maxiter <- control$maxiter
  reltol <- control$reltol
  kappa_method <- control$kappa_method
  nstart <- control$nstart

  # Get kappa solver
  solve_kappa_obj <- movMF:::get_solve_kappa(kappa_method)
  do_kappa <- solve_kappa_obj$do_kappa

  # Set up initialization method
  start_method <- control$start
  if (is.null(start_method)) {
    start_method <- "p"
  }

  # Run EM with multiple random starts
  best_result <- NULL
  best_loglik <- -Inf

  for (start in 1:nstart)
  {
    if (verbose)
      cat(sprintf("\n=== Random start %d/%d ===\n", start, nstart))

    # Initialize parameters
    init <- .initialize_params(x, n, d, G,
                               has_noise, noise_idx,
                               start_method, control)
    mu <- init$mu
    kappa <- init$kappa
    pro <- init$pro

    # EM iterations
    loglik <- .compute_loglik(x, mu, kappa, pro, d, Vinv)
    loglik_old <- loglik - 0.5*abs(loglik)
    iter <- 0
    converged <- FALSE
    degenerate <- FALSE
    loglik_trace <- numeric(maxiter)

    while((loglik - loglik_old)/(1+abs(loglik)) > reltol & iter < maxiter)
    {
      loglik_old <- loglik
      iter <- iter + 1

      # E-step: Compute posterior probabilities
      e_step_result <- .e_step(x, mu, kappa, pro, d, Vinv)
      z <- e_step_result$z

      # M-step: Update parameters
      m_step_result <- try(.m_step(x, z, n, d, do_kappa, has_noise, verbose),
                           silent = TRUE)

      # check for degenerate components
      if (inherits(m_step_result, "try-error"))
      {
        degenerate <- TRUE
        iter <- iter-1
        break
      }
      if (is.null(m_step_result))
      {
        degenerate <- TRUE
        iter <- iter-1
        break
      }

      mu <- m_step_result$mu
      kappa <- m_step_result$kappa
      pro <- m_step_result$pro
      loglik <- .compute_loglik(x, mu, kappa, pro, d, Vinv)
      loglik_trace[iter] <- loglik

      if (verbose)
        cat(sprintf("Iteration %d: loglik = %.*g\n",
                    iter, getOption("digits"), loglik))

    }

    if (iter < maxiter & !degenerate)
    {
      converged <- TRUE
      z <- .e_step(x, mu, kappa, pro, d, Vinv)$z
      if (verbose)
      {
        cat(sprintf("Converged at iteration %d, loglik = %.*g\n",
                    iter, getOption("digits"), loglik))
      }
    }

    # Trim loglik_trace to actual number of iterations
    loglik_trace <- loglik_trace[1:iter]

    # Store result if it's the best so far
    if (!is.na(loglik) && loglik > best_loglik)
    {
      best_loglik <- loglik

      # Transpose mu to be d x G (following mclust convention)
      mu_transposed <- t(mu)

      # Add names to components (1:G or 0:G if noise)
      component_names <- if (has_noise) c(1:G, 0) else 1:G
      component_names_no_noise <- 1:G

      # Name the parameters
      colnames(mu_transposed) <- component_names_no_noise
      names(kappa) <- component_names_no_noise
      names(pro) <- component_names

      # Compute BIC and ICL (following mclust convention: higher is better)
      # Number of parameters:
      # - mu: G * (d-1) (unit vectors have d-1 free parameters)
      # - kappa: G
      # - pro: G (one is determined by sum to 1 constraint, noise prop also counted if present)
      npar <- G * (d - 1) + G + (if (has_noise) G else G - 1)
      bic <- 2 * loglik - npar * log(n)
      C <- matrix(0L, nrow = n, ncol = ncol(z))
      C[cbind(seq_len(n), apply(z, 1, which.max))] = 1
      icl <- bic + 2 * sum(C * ifelse(z > 0, log(z), 0))

      best_result <- list(data = x,
                          n = n,
                          d = d,
                          G = G,
                          parameters = list(mu = mu_transposed,
                                            kappa = kappa,
                                            pro = pro,
                                            Vinv = Vinv),
                          df = npar,
                          # hypvol = if (is.null(Vinv)) as.double(NA) else 1 / Vinv,
                          z = z,
                          classification = map_classification(z, noise = has_noise),
                          loglik = loglik,
                          loglik_trace = loglik_trace,
                          bic = bic,
                          icl = icl,
                          control = control,
                          iterations = iter,
                          converged = converged
      )
    }
  }

  if(verbose)
  {
    cat(sprintf("\n=== Best result ===\n"))
    cat(sprintf("Final log-likelihood: %.*g\n",
                getOption("digits"), best_result$loglik))
  }

  return(best_result)
}

#' Initialize parameters for EM algorithm
#'
#' @param x Data matrix (n x d)
#' @param n Number of observations
#' @param d Number of dimensions
#' @param G Number of components
#' @param has_noise Logical; whether noise component is included
#' @param noise_idx Numeric vector of indices initially assigned to noise (optional)
#' @param start Initialization method: "p" (partition-based), "s"/"S" (seeds), or "i" (random ids)
#' @param control List of control parameters (used to get kappa_method and noise_prop)
#' @return List with mu, kappa, and pro
#' @keywords internal
.initialize_params <- function(x, n, d, G,
                               has_noise = FALSE, noise_idx = NULL,
                               start = "p",
                               control = control_movMFnoise() )
{
  # Get kappa solver from control
  solve_kappa_obj <- movMF:::get_solve_kappa(control$kappa_method)
  do_kappa <- solve_kappa_obj$do_kappa

  # Determine which observations to use for initialization
  if (has_noise && !is.null(noise_idx))
  {
    # Exclude noise observations from initialization
    avail_idx <- setdiff(seq_len(n), noise_idx)
    if (length(avail_idx) < G) {
      stop("Not enough non-noise observations to initialize G components")
    }
    x_init <- x[avail_idx, , drop = FALSE]
    n_init <- length(avail_idx)
  } else {
    x_init <- x
    n_init <- n
  }

  # Initialize posterior probabilities based on method
  if (is.character(start))
  {
    if (start %in% c("p", "S", "s"))
    {
      # Use skmeans-style initialization to get initial centroids
      # Then convert to fuzzy membership with m=2 and cosine dissimilarity
      M <- skmeans:::.skmeans_init_for_normalized_x(x_init, G, start)
      # Compute dissimilarities based on cosine distance
      D <- pmax(1 - tcrossprod(x_init, M), 0)
      # Convert to fuzzy memberships (m=2)
      # Using clue's membership function
      P <- clue:::.memberships_from_cross_dissimilarities(D, 2)
    } else if (start == "i")
    {
      # Initialize by choosing random class ids
      ids <- sample.int(G, n_init, replace = TRUE)
      P <- .posterior_from_ids(ids, G)
    } else {
      stop("Invalid initialization method 'start'")
    }
  } else if (!is.null(dim(start)))
  {
    # Matrix of posteriors provided directly
    P <- start
  } else {
    # Vector of class ids provided
    ids <- match(start, sort(unique(start)))
    P <- .posterior_from_ids(ids, G)
  }

  # Perform one M-step to get initial parameters
  w <- colSums(P)
  pro_init <- w / n_init

  M <- crossprod(P, x_init)
  norms <- sqrt(rowSums(M^2))
  mu_init <- M / ifelse(norms > 0, norms, 1)

  kappa_init <- do_kappa(norms, w, d, nu = 0)

  # Adjust mixing proportions if noise component is present
  if (has_noise)
  {
    # Use noise_prop from control
    noise_prop <- control$noise_prop
    pro_init <- c(pro_init * (1 - noise_prop), noise_prop)
  }

  # Keep mu and kappa separate (rather than folding into theta = kappa * mu)
  # to avoid overflow: for extreme kappa, kappa^2 can overflow to Inf even
  # though kappa and kappa * mu remain finite and well-defined.
  # movMF::movMF follows the same convention, always keeping kappa and the
  # mean directions as separate objects throughout the EM recursions.
  return(list(mu = mu_init, kappa = kappa_init, pro = pro_init))
}

#' Convert class IDs to posterior probability matrix
#'
#' @param ids Vector of class IDs
#' @param k Number of classes
#' @return Binary posterior probability matrix (n x k)
#' @keywords internal
.posterior_from_ids <- function(ids, k)
{
  n <- length(ids)
  P <- matrix(0, nrow = n, ncol = k)
  P[cbind(seq_len(n), ids)] <- 1
  P
}

#' E-step: Compute posterior probabilities
#'
#' @param x Data matrix (n x d)
#' @param mu Mean direction parameters (G x d), each row a unit vector
#' @param kappa Concentration parameters (length G)
#' @param pro Mixing proportions (length G or G+1 if noise)
#' @param d Number of dimensions
#' @param Vinv Inverse hypervolume for noise component (or NULL if no noise)
#' @return List with z (posterior probabilities)
#' @keywords internal
.e_step <- function(x, mu, kappa, pro, d, Vinv = NULL)
{
  G <- nrow(mu)
  n <- nrow(x)
  has_noise <- !is.null(Vinv)

  # Compute cross_prod = x %*% t(kappa * mu) = kappa_g * <x_i, mu_g>
  # NOTE: theta = kappa * mu is formed here (linear in kappa), but kappa
  # itself is never recovered via sqrt(rowSums(theta^2)): squaring kappa
  # can overflow to Inf for extreme concentrations (e.g. kappa ~ 1e200)
  # even though kappa and kappa * mu remain perfectly finite. Always pass
  # kappa in directly (as movMF::movMF does) instead of reconstructing it
  # from theta.
  cross_prod <- tcrossprod(x, kappa * mu)

  # Compute log normalizing constants
  # movMF uses lH(kappa, nu) where nu = d/2-1
  log_C <- -movMF:::lH(kappa, d / 2 - 1)

  # Compute log posteriors: log z[i,g] = log(pro_g) + <x_i, theta_g> + log_C_g
  log_z <- sweep(cross_prod, 2, log(pro[1:G]), "+")
  log_z <- sweep(log_z, 2, log_C, "+")

  # Add noise component if present
  if (has_noise) {
    # Uniform density on sphere: f(x) = Vinv (constant)
    # log f(x) = log(Vinv)
    log_noise <- matrix(log(Vinv) + log(pro[G + 1]), nrow = n, ncol = 1)
    log_z <- cbind(log_z, log_noise)
  }

  # Normalize to get probabilities (softmax)
  max_log_z <- apply(log_z, 1, max)
  z <- exp(log_z - max_log_z)
  z <- z / rowSums(z)

  return(list(z = z))
}

#' M-step: Update parameters
#'
#' @param x Data matrix (n x d)
#' @param z Posterior probabilities (n x G or n x (G+1) if noise)
#' @param n Number of observations
#' @param d Number of dimensions
#' @param do_kappa Function to compute kappa
#' @param has_noise Logical; whether noise component is present
#' @param verbose Logical; print warnings if TRUE
#' @return List with mu, kappa, and pro, or NULL if degenerate
#' @keywords internal
.m_step <- function(x, z, n, d, do_kappa, has_noise = FALSE, verbose)
{
  G <- if (has_noise) ncol(z) - 1 else ncol(z)

  # Extract movMF component posteriors (exclude noise if present)
  z_movMF <- z[, 1:G, drop = FALSE]

  # Update mixing proportions
  w <- colSums(z)
  pro <- w / n

  # Check for degenerate movMF components (very small weight)
  if (any(pro[1:G] < 1e-8)) {
    if (verbose) {
      cat("Warning: Degenerate component detected\n")
    }
    return(NULL)
  }

  # Update mean directions and concentration parameters for each movMF component
  mu <- matrix(0, nrow = G, ncol = d)
  kappa <- numeric(G)

  for (g in 1:G)
  {
    # Compute weighted sum for component g
    M_g <- crossprod(z_movMF[, g], x)

    # Compute norm of weighted sum
    norm_g <- sqrt(sum(M_g^2))

    # Update mu: normalized mean direction
    mu[g, ] <- if (norm_g > 0) M_g / norm_g else M_g

    # Update kappa using the do_kappa function
    kappa[g] <- do_kappa(norm_g, w[g], d, nu = 0)
  }

  return(list(mu = mu, kappa = kappa, pro = pro))
}

#' Compute log-likelihood
#'
#' @param x Data matrix (n x d)
#' @param mu Mean direction parameters (G x d), each row a unit vector
#' @param kappa Concentration parameters (length G)
#' @param pro Mixing proportions (length G or G+1 if noise)
#' @param d Number of dimensions
#' @param Vinv Inverse hypervolume for noise component (or NULL if no noise)
#' @return Log-likelihood value
#' @keywords internal
.compute_loglik <- function(x, mu, kappa, pro, d, Vinv = NULL)
{
  G <- nrow(mu)
  has_noise <- !is.null(Vinv)

  # Compute cross_prod = x %*% t(kappa * mu) = kappa_g * <x_i, mu_g>
  # (see note in .e_step(): kappa is never recovered from squared theta)
  cross_prod <- tcrossprod(x, kappa * mu)

  # Compute log normalizing constants
  # movMF uses lH(kappa, nu) where nu = d/2-1
  log_C <- -movMF:::lH(kappa, d / 2 - 1)

  # Log-likelihood contributions: log(pro_g) + <x_i, theta_g> + log_C_g
  log_densities <- sweep(cross_prod, 2, log(pro[1:G]), "+")
  log_densities <- sweep(log_densities, 2, log_C, "+")

  # Add noise component if present
  if (has_noise) {
    log_noise <- matrix(log(Vinv) + log(pro[G + 1]), nrow = nrow(x), ncol = 1)
    log_densities <- cbind(log_densities, log_noise)
  }

  # Sum over log-sum-exp for each observation
  max_vals <- apply(log_densities, 1, max)
  loglik <- sum(max_vals + log(rowSums(exp(log_densities - max_vals))))

  return(loglik)
}

