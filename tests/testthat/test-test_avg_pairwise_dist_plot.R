data(adults)
data(larvae)

#Test_APD_plot
test_that("avg_pairwise_distance — adults: repeat_number_10 has expected shape", {
  # prevent accidental Rplots.pdf if anything draws
  withr::local_options(device = function(...) grDevices::png(filename = tempfile(fileext = ".png")))
  set.seed(20)
  
  result <- test_threshold(
    adults,
    repeats = 10,
    t_min = 100, t_max = 200, t_step = 10,
    group = "location",
    seed = 20
  )
  
  # sanity: repeat_number_10 exists
  expect_true("repeat_number_10" %in% names(result$avg_distances))
  
  # expected dimensions
  expect_identical(dim(result$avg_distances$repeat_number_10), c(374L, 3L))
})

#Test_APD_plot_threshold_unmet
test_that("avg_pairwise_distance — larvae: warns and repeat_number_10 has expected shape", {
  withr::local_options(device = function(...) grDevices::png(filename = tempfile(fileext = ".png")))
  set.seed(20)
  
  expect_warning({
    result_2 <- test_threshold(
      larvae,
      repeats = 10,
      t_min = 600, t_max = 700, t_step = 10,
      group = "breeding_site_type",
      seed = 20
    )
    
    expect_true("repeat_number_10" %in% names(result_2$avg_distances))
    expect_identical(dim(result_2$avg_distances$repeat_number_10), c(594L, 3L))
  })
})
