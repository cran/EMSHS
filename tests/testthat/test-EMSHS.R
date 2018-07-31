# This file contains tests for assessing the validity of some of the outputs
# produced by EMSHS in the absence and presence of graph information
# (i.e., no edges vs. edges)

context("EMSHS Tests Based on Network Information")
set.seed(100)

y <- matrix(rnorm(5*1), ncol = 1)
X <- matrix(rnorm(5*50), ncol = 50)
mus <- matrix(rnorm(50*1), ncol = 1)
nu <- 1

test_that("EMSHS produces valid outputs in the absence of graph information", {
  E = NULL
  em_no_edge <- EMSHS(y, X, mus, nu, E = NULL,
                      a_sigma = 1, b_sigma = 1, a_omega = 2, b_omega = 1,
                      w = 1, eps = 1e-05)

  # niter should be > 0 and its vector length must be equivalent to the number
  # of shrinkage parameters indicated
  expect_true(all(em_no_edge$niter > 0 & (length(em_no_edge$niter) == length(mus))))

  # beta represents a square matrix of regression coefficients where the rows and cols
  # are equal to the number of predictors
  expect_equal(dim(em_no_edge$beta), c(ncol(X),ncol(X)))

  # in the absence of graph information, there should be no correlation calculations
  # for the shrinkage parameter (omega)
  expect_true(all(length(em_no_edge$omega) == 0))

})

test_that("EMSHS produces valid outputs in the presence of graph information", {
  EE <- matrix(c(26,40,
                 40,26,
                 26,16,
                 16,26,
                 16,1,
                 1,16,
                 1,34,
                 34,1,
                 1,8,
                 8,1,
                 1,21,
                 21,1,
                 34,28,
                 28,34,
                 48,2,
                 2,48,
                 2,3,
                 3,2,
                 3,27,
                 27,3), nrow = 20, ncol = 2, byrow = TRUE)

  E <- EE[do.call(order, lapply(1:ncol(EE), function(i) EE[,i])),]

  em_edge <- EMSHS(y, X, mus, nu, E,
                   a_sigma = 1, b_sigma = 1, a_omega = 2, b_omega = 1,
                   w = 1, eps = 1e-05)

  # niter should be > 0 and its vector length must be equivalent to the number
  # of shrinkage parameters indicated
  expect_true(all(em_edge$niter > 0 & (length(em_edge$niter) == length(mus))))

  # beta represents a square matrix of regression coefficients where the rows and cols
  # are equal to the number of predictors
  expect_equal(dim(em_edge$beta), c(ncol(X),ncol(X)))

  # in the presence of graph info, imputed correlations should be calculated
  # for the shrinkage parameter (omega)
  expect_true(all(length(em_edge$omega) != 0))

})
