test_that("convert.mcmc.list handles coda mcmc and mcmc.list objects", {
  mat1 <- matrix(rnorm(100), nrow = 50, dimnames = list(NULL, c("alpha", "beta")))
  mat2 <- matrix(rnorm(100), nrow = 50, dimnames = list(NULL, c("alpha", "beta")))
  
  m1 <- coda::mcmc(mat1)
  res1 <- convert.mcmc.list(m1)
  expect_s3_class(res1, "mcmc.list")
  expect_equal(length(res1), 1)
  expect_equal(coda::varnames(res1), c("alpha", "beta"))

  ml <- coda::mcmc.list(m1, coda::mcmc(mat2))
  res2 <- convert.mcmc.list(ml)
  expect_identical(res2, ml)
})

test_that("convert.mcmc.list handles plain list of matrices", {
  mat1 <- matrix(rnorm(100), nrow = 50, dimnames = list(NULL, c("theta[1]", "theta[2]")))
  mat2 <- matrix(rnorm(100), nrow = 50, dimnames = list(NULL, c("theta[1]", "theta[2]")))
  raw_list <- list(mat1, mat2)

  res <- convert.mcmc.list(raw_list)
  expect_s3_class(res, "mcmc.list")
  expect_equal(length(res), 2)
  expect_equal(coda::varnames(res), c("theta[1]", "theta[2]"))
})
