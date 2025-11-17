# Tests for movMFnoise package

library(testthat)
library(movMFnoise)

# Test 1: Basic EM algorithm functionality ====================================

test_that("em_movMF fits 2-component mixture without noise", {
  set.seed(123)

  # Generate simple 2-cluster data on 3D sphere
  n <- 200
  x1 <- movMF::rmovMF(100, theta = c(5, 0, 0))
  x2 <- movMF::rmovMF(100, theta = c(0, 5, 0))
  x <- rbind(x1, x2)

  # Fit model
  fit <- em_movMF(x, G = 2, control = control_movMFnoise(nstart = 5))

  # Check structure
  expect_s3_class(fit, "movMFnoise")
  expect_equal(fit$G, 2)
  expect_equal(fit$d, 3)
  expect_equal(fit$n, 200)

  # Check parameters
  expect_equal(length(fit$parameters$pro), 2)
  expect_equal(length(fit$parameters$kappa), 2)
  expect_equal(dim(fit$parameters$mu), c(3, 2))
  expect_null(fit$parameters$Vinv)

  # Check model selection criteria
  expect_true(is.numeric(fit$loglik))
  expect_true(is.numeric(fit$bic))
  expect_true(is.numeric(fit$icl))
  expect_true(fit$bic > 0) # BIC should be positive (mclust convention)

  # Check convergence
  expect_true(is.logical(fit$converged))
  expect_true(fit$iterations > 0)

  # Check classifications
  expect_equal(length(fit$classification), 200)
  expect_true(all(fit$classification %in% 1:2))
})

test_that("em_movMF fits mixture with noise component", {
  set.seed(123)

  # Generate data with concentrated clusters + uniform noise
  n1 <- 100
  n2 <- 100
  n_noise <- 50

  x1 <- movMF::rmovMF(n1, theta = c(10, 0, 0))
  x2 <- movMF::rmovMF(n2, theta = c(0, 10, 0))
  # Uniform noise on sphere
  x_noise <- matrix(rnorm(n_noise * 3), ncol = 3)
  x_noise <- x_noise / sqrt(rowSums(x_noise^2))

  x <- rbind(x1, x2, x_noise)

  # Fit model with noise
  fit <- em_movMF(
    x,
    G = 2,
    noise = TRUE,
    control = control_movMFnoise(nstart = 5)
  )

  # Check structure
  expect_s3_class(fit, "movMFnoise")
  expect_equal(fit$G, 2)
  expect_false(is.null(fit$parameters$Vinv))

  # Check parameters include noise component
  expect_equal(length(fit$parameters$pro), 3) # 2 clusters + noise
  expect_equal(length(fit$parameters$kappa), 2) # Only clusters have kappa
  expect_equal(dim(fit$parameters$mu), c(3, 2)) # Only clusters have mu

  # Check noise component classification (0 = noise in mclust convention)
  expect_true(all(fit$classification %in% 0:2))
  expect_true(sum(fit$classification == 0) > 0) # Some points classified as noise

  # Check posterior probabilities
  expect_equal(dim(fit$z), c(250, 3))
  expect_true(all(abs(rowSums(fit$z) - 1) < 1e-10)) # Row sums = 1
})

# Test 2: fit_movMFnoise automatic model selection ============================

test_that("fit_movMFnoise selects optimal G without noise", {
  set.seed(42)

  # Generate 2-cluster data
  n <- 200
  x1 <- movMF::rmovMF(100, theta = c(5, 0, 0))
  x2 <- movMF::rmovMF(100, theta = c(0, 5, 0))
  x <- rbind(x1, x2)

  # Automatic selection
  fit <- fit_movMFnoise(
    x,
    G = 1:4,
    noise = FALSE,
    control = control_movMFnoise(nstart = 5),
    verbose = FALSE
  )

  # Check structure
  expect_s3_class(fit, "movMFnoise")
  expect_s3_class(fit, "movMFnoise_fitted")

  # Check model selection components
  expect_true(!is.null(fit$best_G))
  expect_true(!is.null(fit$best_noise))
  expect_equal(fit$G, fit$best_G)

  # Check summary table
  expect_true(!is.null(fit$model_summary))
  expect_s3_class(fit$model_summary, "data.frame")
  expect_equal(nrow(fit$model_summary), 4) # 4 values of G
  expect_true(all(
    c("G", "noise", "loglik", "bic", "icl", "converged") %in%
      names(fit$model_summary)
  ))

  # Check all BIC/ICL vectors
  expect_equal(length(fit$bic_all), 4)
  expect_equal(length(fit$icl_all), 4)
  expect_equal(length(fit$loglik_all), 4)

  # Best model should maximize BIC
  expect_equal(fit$bic, max(fit$bic_all, na.rm = TRUE))

  # With 2 clear clusters, should select G=2
  expect_equal(fit$best_G, 2)
})

test_that("fit_movMFnoise compares models with and without noise", {
  set.seed(123)

  # Generate data with true noise component
  n1 <- 80
  n2 <- 80
  n_noise <- 40

  x1 <- movMF::rmovMF(n1, theta = c(15, 0, 0))
  x2 <- movMF::rmovMF(n2, theta = c(0, 15, 0))
  x_noise <- matrix(rnorm(n_noise * 3), ncol = 3)
  x_noise <- x_noise / sqrt(rowSums(x_noise^2))
  x <- rbind(x1, x2, x_noise)

  # Compare with and without noise
  fit <- fit_movMFnoise(
    x,
    G = 1:3,
    noise = c(FALSE, TRUE),
    control = control_movMFnoise(nstart = 5),
    verbose = FALSE
  )

  # Should fit 6 models (3 values of G × 2 noise options)
  expect_equal(nrow(fit$model_summary), 6)

  # Check noise column in summary
  expect_true("noise" %in% names(fit$model_summary))
  expect_true(any(fit$model_summary$noise == TRUE))
  expect_true(any(fit$model_summary$noise == FALSE))

  # Best model should use noise component for this data
  expect_true(fit$best_noise)

  # Check G selection is reasonable (likely G=2 with noise)
  expect_true(fit$best_G %in% 1:3)
})

test_that("fit_movMFnoise uses ICL criterion correctly", {
  set.seed(456)

  # Generate data
  n <- 150
  x1 <- movMF::rmovMF(75, theta = c(8, 0, 0))
  x2 <- movMF::rmovMF(75, theta = c(0, 8, 0))
  x <- rbind(x1, x2)

  # Fit with BIC
  fit_bic <- fit_movMFnoise(
    x,
    G = 1:3,
    criterion = "bic",
    control = control_movMFnoise(nstart = 3),
    verbose = FALSE
  )

  # Fit with ICL
  fit_icl <- fit_movMFnoise(
    x,
    G = 1:3,
    criterion = "icl",
    control = control_movMFnoise(nstart = 3),
    verbose = FALSE
  )

  # Check criterion field
  expect_equal(fit_bic$criterion, "bic")
  expect_equal(fit_icl$criterion, "icl")

  # Both should work and select reasonable models
  expect_true(fit_bic$best_G %in% 1:3)
  expect_true(fit_icl$best_G %in% 1:3)

  # ICL often prefers simpler models (more separated clusters)
  # Both criteria should select based on maximum value
  expect_equal(fit_bic$bic, max(fit_bic$bic_all, na.rm = TRUE))
  expect_equal(fit_icl$icl, max(fit_icl$icl_all, na.rm = TRUE))
})

# Test 3: Methods and utilities ===============================================

test_that("predict.movMFnoise works correctly", {
  set.seed(789)

  # Generate and fit data
  x_train <- rbind(
    movMF::rmovMF(50, theta = c(5, 0, 0)),
    movMF::rmovMF(50, theta = c(0, 5, 0))
  )

  fit <- em_movMF(
    x_train,
    G = 2,
    control = control_movMFnoise(nstart = 3)
  )

  # Generate test data
  x_test <- rbind(
    movMF::rmovMF(10, theta = c(5, 0, 0)),
    movMF::rmovMF(10, theta = c(0, 5, 0))
  )

  # Predict on test data
  pred <- predict(fit, newdata = x_test)

  # Check structure
  expect_type(pred, "list")
  expect_true("z" %in% names(pred))
  expect_true("classification" %in% names(pred))

  # Check dimensions
  expect_equal(nrow(pred$z), 20)
  expect_equal(ncol(pred$z), 2)
  expect_equal(length(pred$classification), 20)

  # Check valid predictions
  expect_true(all(pred$classification %in% 1:2))
  expect_true(all(abs(rowSums(pred$z) - 1) < 1e-10))
})

test_that("summary.movMFnoise works correctly", {
  set.seed(111)

  x <- rbind(
    movMF::rmovMF(50, theta = c(5, 0, 0)),
    movMF::rmovMF(50, theta = c(0, 5, 0))
  )

  fit <- em_movMF(x, G = 2, control = control_movMFnoise(nstart = 3))

  # Get summary
  summ <- summary(fit)

  # Check class
  expect_s3_class(summ, "summary.movMFnoise")

  # Check key components
  expect_true(!is.null(summ$n))
  expect_true(!is.null(summ$d))
  expect_true(!is.null(summ$G))
  expect_true(!is.null(summ$loglik))
  expect_true(!is.null(summ$bic))
  expect_true(!is.null(summ$icl))
  expect_true(!is.null(summ$converged))
  expect_true(!is.null(summ$classification))

  # Should print without error
  expect_output(print(summ), "Mixture of von Mises-Fisher")
})

test_that("print methods work correctly", {
  set.seed(222)

  x <- rbind(
    movMF::rmovMF(50, theta = c(5, 0, 0)),
    movMF::rmovMF(50, theta = c(0, 5, 0))
  )

  # Test em_movMF print
  fit_em <- em_movMF(x, G = 2, control = control_movMFnoise(nstart = 3))
  expect_output(print(fit_em), "movMFnoise")
  expect_output(print(fit_em), "2-component")

  # Test fit_movMFnoise print
  fit_auto <- fit_movMFnoise(
    x,
    G = 1:3,
    control = control_movMFnoise(nstart = 3),
    verbose = FALSE
  )
  expect_output(print(fit_auto), "movMFnoise")
  expect_output(print(fit_auto), "Model selection summary")
  expect_output(print(fit_auto), "Best model")
})

# Test 4: Edge cases and error handling =======================================

test_that("handles single component correctly", {
  set.seed(333)

  x <- movMF::rmovMF(100, theta = c(5, 0, 0))

  fit <- em_movMF(x, G = 1, control = control_movMFnoise(nstart = 3))

  expect_equal(fit$G, 1)
  expect_equal(length(fit$parameters$pro), 1)
  expect_equal(as.numeric(fit$parameters$pro), 1)
  expect_true(all(fit$classification == 1))
})

test_that("handles high-dimensional data", {
  set.seed(444)

  # 10-dimensional sphere
  d <- 10
  n <- 100

  # Generate random direction
  theta1 <- rnorm(d)
  theta1 <- 5 * theta1 / sqrt(sum(theta1^2))
  theta2 <- rnorm(d)
  theta2 <- 5 * theta2 / sqrt(sum(theta2^2))

  x <- rbind(
    movMF::rmovMF(50, theta = theta1),
    movMF::rmovMF(50, theta = theta2)
  )

  fit <- em_movMF(x, G = 2, control = control_movMFnoise(nstart = 3))

  expect_equal(fit$d, 10)
  expect_equal(dim(fit$parameters$mu), c(10, 2))
})

test_that("fit_movMFnoise handles failed models gracefully", {
  set.seed(555)

  # Very small dataset that might cause issues with many components
  x <- movMF::rmovMF(20, theta = c(5, 0, 0))

  # Try to fit many components (some will likely fail or not converge well)
  fit <- fit_movMFnoise(
    x,
    G = 1:5,
    control = control_movMFnoise(nstart = 2, maxiter = 10),
    verbose = FALSE
  )

  # Should still return a result with best model
  expect_s3_class(fit, "movMFnoise_fitted")
  expect_true(!is.null(fit$best_G))
  expect_true(fit$best_G >= 1)
})

# Test 5: BIC/ICL convention (mclust-style) ===================================

test_that("BIC follows mclust convention (positive, maximize)", {
  set.seed(666)

  x <- rbind(
    movMF::rmovMF(50, theta = c(5, 0, 0)),
    movMF::rmovMF(50, theta = c(0, 5, 0))
  )

  # Fit multiple models
  fit1 <- em_movMF(x, G = 1, control = control_movMFnoise(nstart = 3))
  fit2 <- em_movMF(x, G = 2, control = control_movMFnoise(nstart = 3))
  fit3 <- em_movMF(x, G = 3, control = control_movMFnoise(nstart = 3))

  # BIC should be positive (2*loglik - npar*log(n) formula)
  expect_true(fit1$bic > 0)
  expect_true(fit2$bic > 0)
  expect_true(fit3$bic > 0)

  # ICL should also be positive (ICL = BIC + 2*entropy, where entropy can be negative)
  expect_true(fit1$icl > 0)
  expect_true(fit2$icl > 0)
  expect_true(fit3$icl > 0)
  # Note: ICL = BIC + 2*entropy, where entropy <= 0 for well-separated clusters
  # So ICL can be less than or greater than BIC depending on cluster separation

  # Automatic selection should maximize BIC
  fit_auto <- fit_movMFnoise(
    x,
    G = 1:3,
    control = control_movMFnoise(nstart = 3),
    verbose = FALSE
  )

  expect_equal(fit_auto$bic, max(fit_auto$bic_all, na.rm = TRUE))
  best_idx <- which.max(fit_auto$bic_all)
  expect_equal(fit_auto$best_G, fit_auto$G_sequence[best_idx])
})

# Test 6: Model selection validation ==========================================

test_that("fit_movMFnoise correctly identifies true number of components", {
  set.seed(999)

  # Scenario 1: 2 well-separated clusters
  n_per_cluster <- 100
  x1 <- movMF::rmovMF(n_per_cluster, theta = c(15, 0, 0)) # High kappa
  x2 <- movMF::rmovMF(n_per_cluster, theta = c(0, 15, 0)) # High kappa
  x_2clusters <- rbind(x1, x2)

  fit_2 <- fit_movMFnoise(
    x_2clusters,
    G = 1:4,
    control = control_movMFnoise(nstart = 10),
    verbose = FALSE
  )

  # Should identify 2 clusters
  expect_equal(fit_2$best_G, 2)

  # Scenario 2: 3 clusters
  x3 <- movMF::rmovMF(n_per_cluster, theta = c(0, 0, 15))
  x_3clusters <- rbind(x1, x2, x3)

  fit_3 <- fit_movMFnoise(
    x_3clusters,
    G = 1:5,
    control = control_movMFnoise(nstart = 10),
    verbose = FALSE
  )

  # Should identify 3 clusters
  expect_equal(fit_3$best_G, 3)
})

test_that("fit_movMFnoise detects noise component when present", {
  set.seed(777)

  # Data with clear noise
  n_cluster <- 80
  n_noise <- 60

  x1 <- movMF::rmovMF(n_cluster, theta = c(20, 0, 0)) # Very concentrated
  x2 <- movMF::rmovMF(n_cluster, theta = c(0, 20, 0)) # Very concentrated

  # Uniform noise
  x_noise <- matrix(rnorm(n_noise * 3), ncol = 3)
  x_noise <- x_noise / sqrt(rowSums(x_noise^2))

  x_with_noise <- rbind(x1, x2, x_noise)

  # Compare with and without noise
  fit <- fit_movMFnoise(
    x_with_noise,
    G = 1:3,
    noise = c(FALSE, TRUE),
    control = control_movMFnoise(nstart = 10),
    verbose = FALSE
  )

  # Model with noise should be selected
  expect_true(fit$best_noise)

  # Should identify 2 concentrated clusters + noise
  expect_equal(fit$best_G, 2)

  # Check classification recovers structure reasonably well
  # At least some points should be classified as noise
  n_classified_noise <- sum(fit$classification == 0)
  expect_true(n_classified_noise > 0)
  expect_true(n_classified_noise < nrow(x_with_noise))
})
