\name{movMFnoise}
\alias{movMFnoise}
\title{Fit a Mixture of von Mises-Fisher Distributions with Noise}
\description{
Fit a mixture of von Mises-Fisher distributions with a uniform noise component 
on the unit sphere using the EM algorithm.
}
\usage{
movMFnoise(x, G, noise_init = 0.1, max_iter = 100, tol = 1e-6,
           init_method = "kmeans", verbose = FALSE)
}
\arguments{
  \item{x}{Matrix of observations, each row is a point on the unit sphere.}
  \item{G}{Integer. Number of von Mises-Fisher components.}
  \item{noise_init}{Numeric. Initial proportion for the noise component (between 0 and 1).
                    Default is 0.1.}
  \item{max_iter}{Integer. Maximum number of EM iterations. Default is 100.}
  \item{tol}{Numeric. Convergence tolerance for the log-likelihood. Default is 1e-6.}
  \item{init_method}{Character. Initialization method: "kmeans" or "random". 
                     Default is "kmeans".}
  \item{verbose}{Logical. If TRUE, print iteration information. Default is FALSE.}
}
\value{
A list with components:
  \item{theta}{List of estimated parameters (alpha, mu, kappa)}
  \item{noise}{Estimated noise proportion}
  \item{posterior}{Matrix of posterior probabilities (n x (G+1))}
  \item{cluster}{Vector of hard cluster assignments (0 for noise, 1:G for components)}
  \item{loglik}{Vector of log-likelihood values at each iteration}
  \item{iterations}{Number of iterations performed}
  \item{converged}{Logical indicating whether the algorithm converged}
}
\examples{
# Generate data
theta_true <- list(
  alpha = c(0.3, 0.4),
  mu = cbind(c(1, 0, 0), c(0, 1, 0)),
  kappa = c(5, 10)
)
set.seed(123)
sim_data <- rmovMFnoise(n = 300, theta = theta_true)

# Fit model
fit <- movMFnoise(sim_data$data, G = 2)
print(fit$theta)
}
