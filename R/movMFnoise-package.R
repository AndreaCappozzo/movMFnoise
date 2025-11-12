#' movMFnoise: Mixtures of von Mises-Fisher Distributions with Noise
#'
#' @description
#' The movMFnoise package provides functions for fitting mixtures of von Mises-Fisher 
#' distributions with a noise component (uniform density on the sphere). This is useful 
#' for clustering directional data with outliers.
#'
#' @section Main Functions:
#' \describe{
#'   \item{\code{\link{rmovMFnoise}}}{Generate random samples from the mixture model}
#'   \item{\code{\link{dmovMFnoise}}}{Compute density of the mixture model}
#'   \item{\code{\link{movMFnoise}}}{Fit the mixture model using EM algorithm}
#' }
#'
#' @section Model:
#' The model assumes observations come from a mixture of G von Mises-Fisher distributions 
#' plus a uniform noise component:
#' 
#' \deqn{f(x) = \sum_{g=1}^G \alpha_g f_{vMF}(x | \mu_g, \kappa_g) + \alpha_{noise} f_{uniform}(x)}
#' 
#' where \eqn{\alpha_g} are mixing proportions, \eqn{f_{vMF}} is the von Mises-Fisher density,
#' \eqn{\mu_g} are mean directions on the unit sphere, \eqn{\kappa_g} are concentration 
#' parameters, and \eqn{f_{uniform}} is the uniform density on the sphere.
#'
#' @examples
#' # Generate data from a 2-component mixture with noise
#' theta <- list(
#'   alpha = c(0.3, 0.4),
#'   mu = cbind(c(1, 0, 0), c(0, 1, 0)),
#'   kappa = c(5, 10)
#' )
#' set.seed(123)
#' data <- rmovMFnoise(n = 300, theta = theta)
#' 
#' # Fit the model
#' fit <- movMFnoise(data$data, G = 2)
#' print(fit$theta)
#'
#' @docType package
#' @name movMFnoise-package
NULL
