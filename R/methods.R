#' Map observations to components based on posterior probabilities
#'
#' @param z Posterior probability matrix (n x G or n x (G+1) if noise)
#' @param noise Logical indicating if noise component is present (default: FALSE)
#' @return Vector of component assignments (1:G, or 0 for noise if present)
#' @keywords internal
map_classification <- function(z, noise = FALSE)
{
  classification <- apply(z, 1, which.max)

  # If noise component is present, label it as 0 (follows mclust convention)
  # In mclust, noise is labeled as 0
  if (noise) {
    noise_col <- ncol(z)
    # Convert noise assignments from G+1 to 0
    classification[classification == noise_col] <- 0
  }

  return(classification)
}

#' Predict method for movMFnoise objects
#'
#' @param object A `'movMFnoise'` object returned by `movMFnoise()` function
#' call.
#' @param newdata Optional new data matrix (n x d). If missing, training data are
#' used.
#' @param what A character string specifying what to retrieve: `"dens"`
#' returns a vector of values for the mixture density; `"cdens"` returns a
#' matrix of component densities for each mixture component (along the
#' columns); `"z"` returns a matrix of component posterior probabilities; `"map"`
#' returns the map classification.
#' @param logarithm A logical value indicating whether or not the logarithm of
#' the densities/probabilities should be returned.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return Return predictions according to `what` argument.
#' @export
predict.movMFnoise <- function(object, newdata,
                               what = c("dens", "cdens", "z", "map"),
                               normalized = TRUE,
                               logarithm = FALSE, ...)
{
  if(!inherits(object, "movMFnoise"))
    stop("object not of class 'movMFnoise'")
  what <- match.arg(what)
  if(missing(newdata))
  {
    newdata <- object$data
  } else
  {
    newdata <- data.matrix(newdata)
    # Normalize newdata to unit vectors
    newdata <- newdata / sqrt(rowSums(newdata^2))
  }
  n <- nrow(newdata)
  if((d <- ncol(newdata)) != object$d)
    stop("newdata of different dimension from <object>$data")
  G <- object$G
  surface_area <- (2*pi^(d/2)) / gamma(d/2)

  # Extract parameters
  # mu is now d x G, transpose back to G x d for internal use
  mu    <- t(object$parameters$mu)
  kappa <- object$parameters$kappa
  pro   <- object$parameters$pro
  Vinv  <- object$parameters$Vinv
  has_noise <- !is.null(Vinv)

  # Compute cross_prod = x %*% t(kappa * mu) = kappa_g * <x_i, mu_g>
  cross_prod <- tcrossprod(newdata, kappa * mu)
  # Compute log normalizing constants
  log_C <- -movMF:::lH(kappa, d / 2 - 1)
  # Compute log-density = log(pro_g) + <x_i, theta_g> + log_C_g
  ldens <- sweep(cross_prod, 2, log(pro[1:G]), "+")
  ldens <- sweep(ldens, 2, log_C, "+")

  if(has_noise)
  {
    # add noise component if present
    # Uniform density on sphere: f(x) = Vinv (constant)
    # log f(x) = log(Vinv)
    ldens_noise <- matrix(log(Vinv) + log(pro[G + 1]),
                          nrow = n, ncol = 1,
                          dimnames = list(NULL, "0"))
    ldens <- cbind(ldens, ldens_noise)
  }

  if(!normalized)
  {
    ldens <- ldens - log(surface_area)
  }

  if(what == "cdens")
  {
    # return mixture components density
    cdens <- if(logarithm) ldens else exp(ldens)
    return(cdens)
  }

  max_ldens <- apply(ldens, 1, max)
  if(what == "z")
  {
    # return probability of belong to mixture components
    # Normalize to get probabilities (softmax)
    z <- exp(ldens - max_ldens)
    z <- z / rowSums(z)
    return(z)
  }

  if(what == "map")
  {
    # return map classification
    z <- exp(ldens - max_ldens)
    classification <- apply(z, 1, which.max)
    # label noise component as 0 (mclust convention)
    if(has_noise)
      classification[classification == (G + 1)] <- 0
    return(classification)
  }

  # return mixture density
  ldens <- max_ldens + log(rowSums(exp(ldens - max_ldens)))
  dens <- if(logarithm) ldens else exp(ldens)
  return(dens)
}

#' Print method for movMFnoise objects
#'
#' @param x A movMFnoise object
#' @param ... Additional arguments passed to print
#' @export
print.movMFnoise <- function(x, ...)
{
  # Check if this is a fitted model with model selection
  is_fitted <- nrow(x$models_summary) > 1

  txt <- paste0("'movMFnoise' model object: ")
  has_noise <- !is.null(x$parameters$Vinv)

  if (x$G == 0 & has_noise) {
    txt <- paste0(txt, "single noise component")
  } else {
    txt <- paste0(txt, "mixture of ", x$G, " von Mises-Fisher distributions")
    if(has_noise) txt <- paste0(txt, " with noise")
  }
  .catwrap(paste(txt))
  cat("\n")

  # Print model selection summary if available
  if(is_fitted)
  {
    cat("Model selection summary:\n")
    cat(sprintf("  Best model: G=%d, noise=%s (selected by %s)\n",
                x$best_G, x$best_noise, toupper(x$criterion)))
    cat(sprintf("  Number of fitted models: %d\n", nrow(x$models_summary)))
    cat(sprintf("  G range: %d-%d\n",
                min(x$models_summary$G),
                max(x$models_summary$G)))
    cat("\n")
  }

  cat("Available components:\n")
  print(names(x))

  invisible(x)
}

#' Summary method for movMFnoise objects
#'
#' @param object A movMFnoise object
#' @param ... Additional arguments (currently unused)
#' @export
summary.movMFnoise <- function(object, ...)
{
  has_noise <- !is.null(object$parameters$Vinv)
  is_fitted <- nrow(object$models_summary) > 1

  # Get classification (use stored classification if available)
  classification <- if (!is.null(object$classification)) {
    object$classification
  } else {
    # Fallback for older objects without classification component
    map_classification(object$z, noise = has_noise)
  }
  classification <- factor(classification,
                           levels = { l <- seq_len(object$G)
                           if(has_noise) l <- c(l,0)
                           l })

  uncertainty <- 1 - apply(object$z, 1, max)

  # Mixing proportions
  pro <- object$parameters$pro
  if (has_noise) {
    names(pro) <- c(seq_len(object$G), 0)
  } else {
    names(pro) <- seq_len(object$G)
  }

  result <- list(n = object$n,
                 d = object$d,
                 G = object$G,
                 loglik = object$loglik,
                 df = object$df,
                 bic = object$bic,
                 icl = object$icl,
                 pro = pro,
                 mu = object$parameters$mu,
                 kappa = object$parameters$kappa,
                 Vinv = object$parameters$Vinv,
                 hypvol = object$hypvol,
                 classification = classification,
                 uncertainty = uncertainty,
                 converged = object$converged,
                 iterations = object$iterations)

  # Add model selection info if available
  if(is_fitted)
  {
    result$models_summary <- object$models_summary
    result$criterion <- object$criterion
  }

  class(result) <- "summary.movMFnoise"
  return(result)
}

#' Print method for summary.movMFnoise objects
#'
#' @param x A summary.movMFnoise object
#' @param digits Number of digits to print
#' @param ... Additional arguments passed to print
#' @export
print.summary.movMFnoise <- function(x, digits = getOption("digits"), ...)
{
  has_noise <- !is.null(x$Vinv)
  title <- "Mixture of von Mises-Fisher distributions"
  title <- ifelse(has_noise, paste(title, "with noise component"), title)
  txt <- paste(rep("-", min(nchar(x$title), getOption("width"))), collapse = "")
  .catwrap(txt)
  .catwrap(title)
  .catwrap(txt)
  cat("\n")

  tab <- data.frame("n" = x$n, "d" = x$d,
                    "G" = paste0(x$G, if(has_noise) "+noise" else ""),
                    "log-likelihood" = x$loglik, "df" = x$df,
                    "BIC" = x$bic, "ICL" = x$icl,
                    check.names = FALSE)
  print(tab, row.names = FALSE, digits = digits)
  cat("\n")

  cat("Mixing proportions (pi):\n")
  print(round(x$pro, digits = digits))
  cat("\n")

  if (x$G > 0)
  {
    cat("Mean directions (mu):\n")
    .printShortMatrix(x$mu, digits = digits,
                      head = 5, tail = 2, chead = 5, ctail = 2)
    cat("\n")
    cat("Concentration parameters (kappa):\n")
    .printShortVector(x$kappa, digits = digits,
                      head = 5, tail = 2)
    cat("\n")
  }

  if(has_noise)
  {
    .catwrap(sprintf("Noise: uniform distribution on S^%d with density %.*g",
                     x$d-1, getOption("digits"), x$Vinv))
    cat("\n")
  }

  cat("Classification table:")
  print(table(x$classification))

  # Print model comparison table if available
  if(!is.null(x$models_summary))
  {
    cat("\nModel selection:\n")
    summary_table <- x$models_summary
    # mark best model
    best_model <- rep("", nrow(summary_table))
    best_row <- which.max(summary_table[[x$criterion]])
    best_model[best_row] <- "*"
    summary_table <- cbind(summary_table, best_model)
    colnames(summary_table)[ncol(summary_table)] <- ""
    print(summary_table, row.names = FALSE, digits = getOption("digits"))
    cat(sprintf("\bBest model G=%d%s selected by %s.\n",
                x$G, ifelse(is.null(x$Vinv), "", "+noise"),
                toupper(x$criterion)))
  }

  invisible(x)
}
