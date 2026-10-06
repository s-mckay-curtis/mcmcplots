test_that("mcmcplotsPalette returns expected colors and lengths", {
  p_rb1 <- mcmcplotsPalette(1, type = "rainbow")
  expect_equal(length(p_rb1), 1)
  expect_type(p_rb1, "character")

  p_rb4 <- mcmcplotsPalette(4, type = "rainbow")
  expect_equal(length(p_rb4), 4)

  p_seq <- mcmcplotsPalette(3, type = "sequential")
  expect_equal(length(p_seq), 3)

  p_gray <- mcmcplotsPalette(5, type = "grayscale")
  expect_equal(length(p_gray), 5)

  p_cb1 <- mcmcplotsPalette(1, type = "colorblind")
  expect_equal(length(p_cb1), 1)

  p_cb4 <- mcmcplotsPalette(4, type = "colorblind")
  expect_equal(length(p_cb4), 4)

  # Verify default type is colorblind
  expect_equal(mcmcplotsPalette(4), mcmcplotsPalette(4, type = "colorblind"))

  p_vir <- mcmcplotsPalette(6, type = "viridis")
  expect_equal(length(p_vir), 6)
})
