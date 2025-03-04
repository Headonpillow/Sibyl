data(adults)
data(larvae)
test_that("repeated_rarefaction_works", {
  # Test the different combination of graphical variables on the adult dataset
  vdiffr::expect_doppelganger("RR_location_CLOUD_ELLIPSE", repeated_rarefaction(adults, repeats = 10, threshold = 100, group = "location", colorb = "location", cloud = TRUE, ellipse = TRUE, seed=20))
  vdiffr::expect_doppelganger("RR_location_CLOUD", repeated_rarefaction(adults, repeats = 10, threshold = 100, group = "location", colorb = "location", cloud = TRUE, ellipse = FALSE, seed=20))
  vdiffr::expect_doppelganger("RR_location", repeated_rarefaction(adults, repeats = 10, threshold = 100, group = "location", colorb = "location", cloud = FALSE, ellipse = FALSE, seed=20))
  vdiffr::expect_doppelganger("RR_CLOUD_ELLIPSE", repeated_rarefaction(adults, repeats = 10, threshold = 100, group = "sample_id", colorb = "location", cloud = TRUE, ellipse = TRUE, seed=20))
  vdiffr::expect_doppelganger("RR_sample/location_CLOUD", repeated_rarefaction(adults, repeats = 10, threshold = 100, group = "sample_id", colorb = "location", cloud = TRUE, ellipse = FALSE, seed=20))
  vdiffr::expect_doppelganger("RR_sample/location", repeated_rarefaction(adults, repeats = 10, threshold = 100, group = "sample_id", colorb = "location", cloud = FALSE, ellipse = FALSE, seed=20))
  #  Coloring the ellipses by group is in the following case impossible. Could be implemented in the future.
  #  vdiffr::expect_doppelganger("RR_sample/sample_CLOUD_ELLIPSE", repeated_rarefaction(adults, repeats = 10, threshold = 100, group = "location", colorb = "sample_id", cloud = TRUE, ellipse = TRUE))
  vdiffr::expect_doppelganger("RR_sample/sample_CLOUD", repeated_rarefaction(adults, repeats = 10, threshold = 100, group = "location", colorb = "sample_id", cloud = TRUE, ellipse = FALSE, seed=20))
  vdiffr::expect_doppelganger("RR_sample/sample", repeated_rarefaction(adults, repeats = 10, threshold = 100, group = "location", colorb = "sample_id", cloud = FALSE, ellipse = FALSE, seed=20))
  # Test behavior when rarefaction threshold is not met on the larval dataset
  testthat::expect_warning(vdiffr::expect_doppelganger("RR_threshold_unmet", repeated_rarefaction(larvae, repeats = 10, threshold = 700, group = "breeding_site_type", colorb = "breeding_site_type", cloud = TRUE, ellipse = TRUE, seed=20)))
})
