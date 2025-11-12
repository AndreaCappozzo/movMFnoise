\name{dmovMFnoise}
\alias{dmovMFnoise}
\title{Density of a Mixture of von Mises-Fisher Distributions with Noise}
\description{
Compute the density of a mixture of von Mises-Fisher distributions with 
a uniform noise component on the unit sphere.
}
\usage{
dmovMFnoise(x, theta, noise = NULL, log = FALSE)
}
\arguments{
  \item{x}{Matrix of observations, each row is a point on the unit sphere.}
  \item{theta}{List containing the mixture model parameters:
    \describe{
      \item{\code{alpha}}{Numeric vector of mixing proportions for the von Mises-Fisher components}
      \item{\code{mu}}{Matrix where each column is a mean direction on the unit sphere 
            for each component}
      \item{\code{kappa}}{Numeric vector of concentration parameters for each component}
    }
  }
  \item{noise}{Numeric. Proportion of the uniform noise component. 
               Default is \code{1 - sum(alpha)}.}
  \item{log}{Logical. If TRUE, return log-density. Default is FALSE.}
}
\value{
Numeric vector of density values for each observation.
}
\examples{
# Compute density for some points
theta <- list(
  alpha = c(0.4, 0.4),
  mu = cbind(c(1, 0, 0), c(0, 1, 0)),
  kappa = c(5, 10)
)
x <- rbind(c(1, 0, 0), c(0, 1, 0), c(0, 0, 1))
dmovMFnoise(x, theta)
}
