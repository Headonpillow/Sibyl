data(adults)
data(larvae)
test_that("APD_calculation_works", {
  # Test normal behavior
  result <- test_threshold(adults, repeats = 10, t_min = 100, t_max = 200, t_step = 10, group = "location", seed=20)
  vdiffr::expect_doppelganger("Test_APD_plot", avg_pairwise_dist_plot(result$avg_distances$repeat_number_10))
  # Test behavior when a threshold is not met
  testthat::expect_warning(result_2 <- test_threshold(larvae, repeats = 10, t_min = 600, t_max = 700, t_step = 10, group = "breeding_site_type", seed=20))
  vdiffr::expect_doppelganger("Test_APD_plot_threshold_unmet", avg_pairwise_dist_plot(result_2$avg_distances$repeat_number_10))
})