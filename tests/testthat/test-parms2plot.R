test_that("parms2plot matches parameter names correctly", {
  varnames <- c("alpha[1]", "alpha[2]", "beta[1]", "beta[2]", "sigma", "gamma")

  # Match all by default (returns all, grouped by scalar vs leaf)
  all_p <- parms2plot(varnames)
  expect_setequal(unname(unlist(all_p)), varnames)

  # Match exact parameter group
  alpha_p <- parms2plot(varnames, parms = "alpha")
  expect_equal(unname(unlist(alpha_p)), c("alpha[1]", "alpha[2]"))

  # Match regex
  reg_p <- parms2plot(varnames, regex = "^(beta|sigma)")
  expect_equal(unname(unlist(reg_p)), c("beta[1]", "beta[2]", "sigma"))

  # Random subset
  set.seed(42)
  rand_p <- parms2plot(varnames, parms = "alpha", random = 1)
  expect_equal(length(unlist(rand_p)), 1)
  expect_true(unlist(rand_p) %in% c("alpha[1]", "alpha[2]"))
})
