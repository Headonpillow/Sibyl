data(adults)

test_that("Accumulation_works (title only)", {
  # prevent Rplots.pdf during non-interactive tests
  withr::local_options(device = function(...) grDevices::png(filename = tempfile(fileext = ".png")))
  
  p <- accumulation_test(adults, step = 50)
  
  # Sanity: structure is as expected
  expect_true(!is.null(p$accumulation_plot), "accumulation_plot is missing")
  expect_s3_class(p$accumulation_plot, "ggplot")
  
  # Title check
  expect_identical(p$accumulation_plot$labels$title, "Species Accumulation Curves")
})
