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
#' @param maxiter Integer. Maximum number of EM iterations. Default is 1000.
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
#' @param noise_prop Proportion of noise used for initialization of EM algorithm.
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
  maxiter = 1000L,
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

#' @keywords internal
#' @noRd
.catwrap <- function(x, width = getOption("width"), ...)
{
# version of cat with wrapping at specified width
  cat(paste(strwrap(x, width = width, ...), collapse = "\n"), "\n")
}

#' Print Representative Columns and Rows for a Matrix
#'
#' Prints a shortened matrix containing configurable numbers of leading and
#' trailing rows and columns.
#'
#' @param x A matrix to summarize.
#' @param head number of initial rows to print
#' @param tail number of last rows to print
#' @param chead number of initial columns to print
#' @param ctail number of last columns to print
#' @param ... extra arguments passed to [print()]
#' @keywords internal
#' @noRd
.printShortMatrix <- function(x, head = 2, tail = 1, chead = 5, ctail = 1, ...)
{
  x <- as.matrix(x)
  nr <- nrow(x)
  nc <- ncol(x)
  rnames <- rownames(x)
  cnames <- colnames(x)
  dnames <- names(dimnames(x))

  if(is.na(head <- as.numeric(head))) head <- 2
  if(is.na(tail <- as.numeric(tail))) tail <- 1
  if(is.na(chead <- as.numeric(chead))) chead <- 5
  if(is.na(ctail <- as.numeric(ctail))) ctail <- 1

  if(nr > (head + tail))
  {
    if(is.null(rnames))
      rnames <- paste("[", 1:nr, ",]", sep ="")
    x <- rbind(x[1:head,,drop=FALSE],
               rep(NA, nc),
               x[(nr-tail+1):nr,,drop=FALSE])
    rownames(x) <- c(rnames[1:head], ":", rnames[(nr-tail+1):nr])
  }
  if(nc > (chead + ctail))
  {
    if(is.null(cnames))
      cnames <- paste("[,", 1:nc, "]", sep ="")
    x <- cbind(x[,1:chead,drop=FALSE],
               rep(NA, nrow(x)),
               x[,(nc-ctail+1):nc,drop=FALSE])
    colnames(x) <- c(cnames[1:chead], "...", cnames[(nc-ctail+1):nc])
  }
  names(dimnames(x)) <- dnames
  print(x, na.print = "", ...)
  invisible(x)
}

#' Print representative values of a vector
#'
#' Prints a shortened vector containing configurable numbers of leading and
#' trailing values
#'
#' @param x A vector to show.
#' @param head number of initial values to print.
#' @param tail number of last values to print.
#' @param ... extra arguments passed to [print()]
#' @keywords internal
#' @noRd
.printShortVector <- function(x, head = 2, tail = 1, ...)
{
  names <- names(x)
  x <- setNames(as.vector(x), names)
  n <- length(x)

  if(is.na(head <- as.numeric(head))) head <- 2
  if(is.na(tail <- as.numeric(tail))) tail <- 1

  if(n > (head + tail))
  {
    if(is.null(names))
      names <- paste("[", 1:n, "]", sep ="")
    x <- c(x[1:head, drop=FALSE], NA,
           x[(n-tail+1):n, drop=FALSE])
    names(x) <- c(names[1:head], "...", names[(n-tail+1):n])
  }
  print.default(x, na.print = "", ...)
  invisible(x)
}
