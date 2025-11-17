#' Map observations to components based on posterior probabilities
#'
#' @param z Posterior probability matrix (n x G or n x (G+1) if noise)
#' @param noise Logical indicating if noise component is present (default: FALSE)
#' @return Vector of component assignments (1:G, or 0 for noise if present)
#' @keywords internal
map_classification <- function(z, noise = FALSE) {
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
#' @param object A movMFnoise object returned by em_movMF
#' @param newdata Optional new data matrix (n x d). If missing, uses training data
#' @param ... Additional arguments (currently unused)
#' @return List with components:
#'   \item{classification}{Vector of component assignments}
#'   \item{z}{Posterior probability matrix}
#' @export
predict.movMFnoise <- function(object, newdata, ...) {
  if (missing(newdata)) {
    # Use the posterior probabilities from fitting
    z <- object$z
  } else {
    # Normalize newdata to unit vectors
    newdata <- newdata / sqrt(rowSums(newdata^2))

    # Extract parameters
    # mu is now d x G, transpose back to G x d for internal use
    mu_internal <- t(object$parameters$mu)
    theta <- mu_internal * object$parameters$kappa
    pro <- object$parameters$pro
    d <- object$d
    Vinv <- object$parameters$Vinv

    # Compute posterior probabilities for new data
    z <- .e_step(newdata, theta, pro, d, Vinv)$z
  }

  # Map to hard classifications
  has_noise <- !is.null(object$parameters$Vinv)
  G <- object$G

  classification <- apply(z, 1, which.max)

  # Label noise component as 0 (mclust convention)
  if (has_noise) {
    classification[classification == (G + 1)] <- 0
  }

  # Set column names for z matrix
  colnames(z) <- if (has_noise) c(seq_len(G), 0) else seq_len(G)

  return(list(classification = classification, z = z))
}

#' Print method for movMFnoise objects
#'
#' @param x A movMFnoise object
#' @param ... Additional arguments passed to print
#' @export
print.movMFnoise <- function(x, ...) {
  # Check if this is a fitted model with model selection
  is_fitted <- inherits(x, "movMFnoise_fitted")

  txt <- paste0("'movMFnoise' model object: ")
  has_noise <- !is.null(x$parameters$Vinv)

  if (x$G == 0 & has_noise) {
    txt <- paste0(txt, "single noise component")
  } else {
    txt <- paste0(
      txt,
      x$G,
      "-component mixture of von Mises-Fisher distributions"
    )
    if (has_noise) {
      txt <- paste0(txt, " with noise")
    }
  }

  cat(txt, "\n")

  # Print model selection summary if available
  if (is_fitted && !is.null(x$model_summary)) {
    cat("\nModel selection summary:\n")
    cat(sprintf(
      "  Best model: G=%d, noise=%s (selected by %s)\n",
      x$best_G,
      x$best_noise,
      toupper(x$criterion)
    ))
    cat(sprintf("  Number of models fitted: %d\n", nrow(x$model_summary)))
    cat(sprintf("  G range: %d-%d\n", min(x$G_sequence), max(x$G_sequence)))
  }

  cat("\nAvailable components:\n")
  print(names(x))

  invisible(x)
}

#' Summary method for movMFnoise objects
#'
#' @param object A movMFnoise object
#' @param ... Additional arguments (currently unused)
#' @export
summary.movMFnoise <- function(object, ...) {
  has_noise <- !is.null(object$parameters$Vinv)
  is_fitted <- inherits(object, "movMFnoise_fitted")

  # Get classification (use stored classification if available)
  classification <- if (!is.null(object$classification)) {
    object$classification
  } else {
    # Fallback for older objects without classification component
    map_classification(object$z, noise = has_noise)
  }
  uncertainty <- 1 - apply(object$z, 1, max)

  # Mixing proportions
  pro <- object$parameters$pro
  if (has_noise) {
    names(pro) <- c(seq_len(object$G), 0)
  } else {
    names(pro) <- seq_len(object$G)
  }

  result <- list(
    n = object$n,
    d = object$d,
    G = object$G,
    loglik = object$loglik,
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
    iterations = object$iterations
  )

  # Add model selection info if available
  if (is_fitted) {
    result$model_summary <- object$model_summary
    result$best_G <- object$best_G
    result$best_noise <- object$best_noise
    result$criterion <- object$criterion
  }

  class(result) <- c(
    "summary.movMFnoise",
    if (is_fitted) "summary.movMFnoise_fitted"
  )

  return(result)
}

#' Print method for summary.movMFnoise objects
#'
#' @param x A summary.movMFnoise object
#' @param digits Number of digits to print
#' @param ... Additional arguments passed to print
#' @export
print.summary.movMFnoise <- function(x, digits = getOption("digits"), ...) {
  has_noise <- !is.null(x$Vinv)
  is_fitted <- inherits(x, "summary.movMFnoise_fitted")

  cat(paste(rep("-", 70), collapse = ""), "\n")
  cat("Mixture of von Mises-Fisher distributions")
  if (has_noise) {
    cat(" with noise component")
  }
  cat("\n")
  cat(paste(rep("-", 70), collapse = ""), "\n\n")

  cat(sprintf("Number of observations: %d\n", x$n))
  cat(sprintf("Dimension: %d\n", x$d))
  cat(sprintf("Number of components: %d", x$G))
  if (has_noise) {
    cat(" + noise")
  }
  cat("\n")

  # Print model selection info if available
  if (is_fitted && !is.null(x$model_summary)) {
    cat(sprintf(
      "Model selection: Best G=%d (noise=%s) selected by %s\n",
      x$best_G,
      x$best_noise,
      toupper(x$criterion)
    ))
  }

  cat(sprintf("Log-likelihood: %.4f\n", x$loglik))
  cat(sprintf("BIC: %.4f\n", x$bic))
  cat(sprintf("ICL: %.4f\n", x$icl))
  cat(sprintf("Converged: %s (iterations: %d)\n\n", x$converged, x$iterations))

  cat("Mixing proportions:\n")
  print(round(x$pro, digits = digits))
  cat("\n")

  if (x$G > 0) {
    cat("Concentration parameters (kappa):\n")
    print(round(x$kappa, digits = digits))
    cat("\n")
  }

  if (has_noise) {
    cat(sprintf("Noise component hypervolume: %.4e\n", x$hypvol))
    cat(sprintf("Noise component density: %.4e\n\n", x$Vinv))
  }

  cat("Classification table:\n")
  print(table(x$classification))
  cat("\n")

  # Print model comparison table if available
  if (is_fitted && !is.null(x$model_summary)) {
    cat("All fitted models:\n")
    summary_table <- x$model_summary
    summary_table$loglik <- round(summary_table$loglik, 2)
    summary_table$bic <- round(summary_table$bic, 2)
    summary_table$icl <- round(summary_table$icl, 2)

    # Mark best model
    best_row <- which(
      summary_table$G == x$best_G &
        summary_table$noise == x$best_noise
    )

    print(summary_table, row.names = FALSE)
    cat(sprintf("  (Row %d selected as best)\n\n", best_row))
  }

  invisible(x)
}
