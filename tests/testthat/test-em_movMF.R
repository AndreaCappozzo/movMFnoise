library(testthat)

test_that("em_movMF basic functionality", {
  skip_if_not_installed("movMF")
  
  # Create simple test data on unit sphere
  set.seed(123)
  n <- 100
  d <- 3
  
  # Generate data from two vMF distributions
  # Component 1: centered at (1, 0, 0)
  mu1 <- c(1, 0, 0)
  # Component 2: centered at (0, 1, 0)
  mu2 <- c(0, 1, 0)
  
  # Simple synthetic data (not perfectly from vMF, but on unit sphere)
  n1 <- 50
  n2 <- 50
  
  # Generate points near mu1
  x1 <- matrix(rnorm(n1 * d), n1, d)
  x1[, 1] <- x1[, 1] + 3  # Bias towards mu1
  x1 <- x1 / sqrt(rowSums(x1^2))  # Normalize to unit sphere
  
  # Generate points near mu2
  x2 <- matrix(rnorm(n2 * d), n2, d)
  x2[, 2] <- x2[, 2] + 3  # Bias towards mu2
  x2 <- x2 / sqrt(rowSums(x2^2))  # Normalize to unit sphere
  
  x <- rbind(x1, x2)
  
  # Test with G = 2
  result <- em_movMF(x, G = 2, control = list(nstart = 1, maxiter = 50, verbose = FALSE))
  
  # Check structure
  expect_true(is.list(result))
  expect_true("parameters" %in% names(result))
  expect_true("z" %in% names(result))
  expect_true("loglik" %in% names(result))
  expect_true("n" %in% names(result))
  expect_true("d" %in% names(result))
  expect_true("G" %in% names(result))
  
  # Check dimensions
  expect_equal(result$n, n)
  expect_equal(result$d, d)
  expect_equal(result$G, 2)
  
  # Check parameter dimensions
  expect_equal(nrow(result$parameters$mu), 2)
  expect_equal(ncol(result$parameters$mu), d)
  expect_equal(length(result$parameters$kappa), 2)
  expect_equal(length(result$parameters$pro), 2)
  
  # Check z dimensions
  expect_equal(nrow(result$z), n)
  expect_equal(ncol(result$z), 2)
  
  # Check that mixing proportions sum to 1
  expect_equal(sum(result$parameters$pro), 1, tolerance = 1e-10)
  
  # Check that z rows sum to 1
  expect_true(all(abs(rowSums(result$z) - 1) < 1e-10))
  
  # Check that mu vectors are unit vectors
  mu_norms <- sqrt(rowSums(result$parameters$mu^2))
  expect_true(all(abs(mu_norms - 1) < 1e-10))
  
  # Check that kappa values are positive
  expect_true(all(result$parameters$kappa > 0))
  
  # Check that class is set
  expect_true(inherits(result, "movMF"))
})

test_that("em_movMF handles control parameters", {
  skip_if_not_installed("movMF")
  
  set.seed(456)
  n <- 50
  d <- 3
  x <- matrix(rnorm(n * d), n, d)
  x <- x / sqrt(rowSums(x^2))
  
  # Test with custom control parameters
  result <- em_movMF(x, G = 1, control = list(
    maxiter = 20,
    nstart = 2,
    verbose = FALSE
  ))
  
  expect_equal(result$control$maxiter, 20)
  expect_equal(result$control$nstart, 2)
  expect_false(result$control$verbose)
})

test_that("em_movMF handles single component", {
  skip_if_not_installed("movMF")
  
  set.seed(789)
  n <- 30
  d <- 3
  x <- matrix(rnorm(n * d), n, d)
  x <- x / sqrt(rowSums(x^2))
  
  result <- em_movMF(x, G = 1, control = list(nstart = 1, verbose = FALSE))
  
  expect_equal(result$G, 1)
  expect_equal(nrow(result$parameters$mu), 1)
  expect_equal(length(result$parameters$kappa), 1)
  expect_equal(length(result$parameters$pro), 1)
  expect_equal(result$parameters$pro[1], 1)
})
