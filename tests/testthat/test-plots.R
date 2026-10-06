# Helper to generate reproducible fake mcmc output
make_test_mcmc <- function() {
  nc <- 3
  nr <- 100
  pnames <- c("alpha[1]", "alpha[2]", "beta[1]", "beta[2]")
  coda::as.mcmc.list(
    lapply(1:nc, function(i) {
      coda::mcmc(matrix(rnorm(nr * length(pnames)), nrow = nr, dimnames = list(NULL, pnames)))
    })
  )
}

test_that("caterplot runs without error", {
  m <- make_test_mcmc()
  pdf(NULL)
  on.exit(dev.off())
  expect_no_error(caterplot(m, "alpha"))
  expect_no_error(caterplot(m, collapse = FALSE))
})

test_that("denplot and traplot run without error", {
  m <- make_test_mcmc()
  pdf(NULL)
  on.exit(dev.off())
  expect_no_error(denplot(m, "alpha", style = "plain"))
  expect_no_error(denplot(m, "alpha", style = "gray"))
  expect_no_error(traplot(m, "alpha", style = "plain"))
  expect_no_error(traplot(m, "alpha", style = "gray"))
})

test_that("mcmcplot generates HTML report and image files", {
  m <- make_test_mcmc()
  td <- tempfile("mcmcplot_test")
  dir.create(td)
  on.exit(unlink(td, recursive = TRUE))

  res <- mcmcplot(m, parms = "alpha", dir = td, filename = "test_report", browse = FALSE)
  html_file <- file.path(td, "test_report.html")
  expect_true(file.exists(html_file))

  # Check that PNG files were produced in directory
  png_files <- list.files(td, pattern = "\\.png$")
  expect_gt(length(png_files), 0)

  # Check that HTML file contains title and images
  html_content <- readLines(html_file)
  expect_true(any(grepl("MCMC Plots", html_content)))
  expect_true(any(grepl("<img", html_content)))
})
