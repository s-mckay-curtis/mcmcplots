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

test_that("caterplot runs without error across styles", {
  m <- make_test_mcmc()
  pdf(NULL)
  on.exit(dev.off())
  expect_no_error(caterplot(m, "alpha", style = "clean"))
  expect_no_error(caterplot(m, "alpha", style = "plain"))
  expect_no_error(caterplot(m, "alpha", style = "gray"))
  expect_no_error(caterplot(m, collapse = FALSE, style = "clean"))
})

test_that("denplot, traplot, and rmeanplot run without error across styles", {
  m <- make_test_mcmc()
  pdf(NULL)
  on.exit(dev.off())
  expect_no_error(denplot(m, "alpha", style = "clean"))
  expect_no_error(denplot(m, "alpha", style = "plain"))
  expect_no_error(denplot(m, "alpha", style = "gray"))
  expect_no_error(traplot(m, "alpha", style = "clean"))
  expect_no_error(traplot(m, "alpha", style = "plain"))
  expect_no_error(traplot(m, "alpha", style = "gray"))
  expect_no_error(rmeanplot(m, "alpha", style = "clean"))
  expect_no_error(rmeanplot(m, "alpha", style = "plain"))
  expect_no_error(rmeanplot(m, "alpha", style = "gray"))
})

test_that("mcmcplot1 renders with clean style and custom colors", {
  m <- make_test_mcmc()
  pdf(NULL)
  on.exit(dev.off())
  expect_no_error(mcmcplot1(m[, "alpha[1]", drop = FALSE], style = "clean"))
  expect_no_error(mcmcplot1(m[, "alpha[1]", drop = FALSE], style = "plain"))
  expect_no_error(mcmcplot1(m[, "alpha[1]", drop = FALSE], style = "gray"))
  # Verify custom col works and flows through to rmeanplot1 without error
  custom_cols <- c("firebrick", "darkblue", "goldenrod")
  expect_no_error(mcmcplot1(m[, "alpha[1]", drop = FALSE], col = custom_cols))
})

test_that("mcmcplot generates modern HTML report with search and cards", {
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

  # Check modern HTML structure
  html_content <- readLines(html_file)
  expect_true(any(grepl("<!DOCTYPE html>", html_content)))
  expect_true(any(grepl("id=\"param_search\"", html_content)))
  expect_true(any(grepl("class=\"plot-card\"", html_content)))
  expect_true(any(grepl("function filterPlots", html_content)))
  expect_true(any(grepl("<img", html_content)))
})

test_that("mcmcplot supports retina scaling and base64 self-contained embedding", {
  m <- make_test_mcmc()
  td <- tempfile("mcmcplot_test_b64")
  dir.create(td)
  on.exit(unlink(td, recursive = TRUE))

  res <- mcmcplot(m, parms = "alpha[1]", dir = td, filename = "b64_report",
                  browse = FALSE, retina = TRUE, embed.img = TRUE)
  html_file <- file.path(td, "b64_report.html")
  expect_true(file.exists(html_file))

  html_content <- paste(readLines(html_file), collapse = "\n")
  expect_true(grepl("data:image/png;base64,", html_content))

  # Test custom pointsize
  expect_no_error(mcmcplot(m, parms = "alpha[1]", dir = td, filename = "ps_report",
                           browse = FALSE, pointsize = 16))
})
