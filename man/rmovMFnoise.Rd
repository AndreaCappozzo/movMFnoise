\name{rmovMFnoise}
\alias{rmovMFnoise}
\title{Generate Random Samples from a Mixture of von Mises-Fisher Distributions with Noise}
\description{
Generate random samples from a mixture of von Mises-Fisher distributions with 
a uniform noise component on the unit sphere.
}
\usage{
rmovMFnoise(n, theta, noise = NULL)
}
\arguments{
  \item{n}{Integer. Number of observations to generate.}
  \item{theta}{List containing the mixture model parameters:
    \describe{
      \item{\code{alpha}}{Numeric vector of mixing proportions for the von Mises-Fisher 
            components (should sum to less than or equal to 1)}
      \item{\code{mu}}{Matrix where each column is a mean direction on the unit sphere 
            for each component}
      \item{\code{kappa}}{Numeric vector of concentration parameters (>= 0) for each component}
    }
  }
  \item{noise}{Numeric. Proportion of observations from the uniform noise component 
               (between 0 and 1). Default is 0, meaning \code{1 - sum(alpha)}.}
}
\value{
A list with components:
  \item{data}{Matrix of generated observations, each row is a point on the unit sphere}
  \item{cluster}{Integer vector indicating the true cluster membership 
         (0 for noise, 1:G for components)}
}
\examples{
# Generate data from 2 components plus noise
theta <- list(
  alpha = c(0.3, 0.4),
  mu = cbind(c(1, 0, 0), c(0, 1, 0)),
  kappa = c(5, 10)
)
set.seed(123)
data <- rmovMFnoise(n = 300, theta = theta)
table(data$cluster)
}
