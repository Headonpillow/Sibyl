withr::local_seed(1234)
data(adults)
data(larvae)
test_that("repeated_rarefaction_works", {
  # Test normal behavior with a single threshold
  result <- test_threshold(adults, repeats = 10, t_min = 100, t_max = 200, t_step = 10, group = "location")
  vdiffr::expect_doppelganger("TT_single_threshold", result$index_plot)
  # Test normal behavior with multiple repeat numbers
  result_2 <- test_threshold(adults, repeats = c(10,100), t_min = 100, t_max = 200, t_step = 10, group = "location")
  vdiffr::expect_doppelganger("TT_multiple_threshold", result_2$index_plot)
  # Test behavior when samples do not meet the threshold (for index and ordinations)
  result_3 <- test_threshold(larvae, repeats = 10, t_min = 600, t_max = 700, t_step = 10, group = "breeding_site_type")
  testthat::expect_warning(vdiffr::expect_doppelganger("TT_multiple_threshold_unmet_index", result_3$index_plot))
  testthat::expect_warning(vdiffr::expect_doppelganger("TT_multiple_threshold_unmet_ordination", result_3$ordination_plots$`repeat_number 10`[[1]]))
})
