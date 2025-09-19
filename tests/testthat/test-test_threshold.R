data(adults)
data(larvae)

# Test normal behavior with a single threshold
test_that("test_threshold — adults: single repeat index_plot shape", {
  withr::local_options(device = function(...) grDevices::png(filename = tempfile(fileext = ".png")))
  set.seed(20)
  
  result <- test_threshold(
    adults,
    repeats = 10,
    t_min = 100, t_max = 200, t_step = 10,
    group = "location",
    seed = 20
  )
  
  expect_true(!is.null(result$index_plot$data))
  expect_identical(dim(result$index_plot$data), c(11L, 3L))
})

# Test normal behavior with multiple repeat numbers
test_that("test_threshold — adults: multiple thresholds index_plot shape", {
  withr::local_options(device = function(...) grDevices::png(filename = tempfile(fileext = ".png")))
  set.seed(20)
  
  result_2 <- test_threshold(
    adults,
    repeats = c(10, 100),
    t_min = 100, t_max = 200, t_step = 10,
    group = "location",
    seed = 20
  )
  
  expect_true(!is.null(result_2$index_plot$data))
  expect_identical(dim(result_2$index_plot$data), c(22L, 3L))
})

# Test behavior when samples do not meet the threshold (for index plot)
test_that("test_threshold — larvae: threshold unmet warns and index_plot shape", {
  withr::local_options(device = function(...) grDevices::png(filename = tempfile(fileext = ".png")))
  set.seed(20)
  
  expect_warning({
    result_3 <- test_threshold(
      larvae,
      repeats = 10,
      t_min = 600, t_max = 700, t_step = 10,
      group = "breeding_site_type",
      seed = 20
    )
    
    expect_true(!is.null(result_3$index_plot$data))
    expect_identical(dim(result_3$index_plot$data), c(11L, 3L))
  })
})
