withr::local_seed(1234)
data(adults)
data(larvae)
test_that("repeated_rarefaction_works", {
  # Test the different combination of graphical variables on the adult dataset
  vdiffr::expect_doppelganger("Repeated_rarefaction_location_CLOUD_ELLIPSE", repeated_rarefaction(adults, repeats = 10, threshold = 100, group = "location", colorb = "location", cloud = TRUE, ellipse = TRUE))
  vdiffr::expect_doppelganger("Repeated_rarefaction_location_CLOUD", repeated_rarefaction(adults, repeats = 10, threshold = 100, group = "location", colorb = "location", cloud = TRUE, ellipse = FALSE))
  vdiffr::expect_doppelganger("Repeated_rarefaction_location", repeated_rarefaction(adults, repeats = 10, threshold = 100, group = "location", colorb = "location", cloud = FALSE, ellipse = FALSE))
  vdiffr::expect_doppelganger("Repeated_rarefaction_sample/location_CLOUD_ELLIPSE", repeated_rarefaction(adults, repeats = 10, threshold = 100, group = "sample_id", colorb = "location", cloud = TRUE, ellipse = TRUE))
  vdiffr::expect_doppelganger("Repeated_rarefaction_sample/location_CLOUD", repeated_rarefaction(adults, repeats = 10, threshold = 100, group = "sample_id", colorb = "location", cloud = TRUE, ellipse = FALSE))
  vdiffr::expect_doppelganger("Repeated_rarefaction_sample/location", repeated_rarefaction(adults, repeats = 10, threshold = 100, group = "sample_id", colorb = "location", cloud = FALSE, ellipse = FALSE))
  vdiffr::expect_doppelganger("Repeated_rarefaction_sample/sample_CLOUD_ELLIPSE", repeated_rarefaction(adults, repeats = 10, threshold = 100, group = "location", colorb = "sample_id", cloud = TRUE, ellipse = TRUE))
  vdiffr::expect_doppelganger("Repeated_rarefaction_sample/sample_CLOUD", repeated_rarefaction(adults, repeats = 10, threshold = 100, group = "location", colorb = "sample_id", cloud = TRUE, ellipse = FALSE))
  vdiffr::expect_doppelganger("Repeated_rarefaction_sample/sample", repeated_rarefaction(adults, repeats = 10, threshold = 100, group = "location", colorb = "sample_id", cloud = FALSE, ellipse = FALSE))
  # Test behavior when rarefaction threshold is not met on the larval dataset
  vdiffr::expect_doppelganger("Repeated_rarefaction_threshold_unmet", repeated_rarefaction(larvae, repeats = 10, threshold = 700, group = "breeding_site_type", colorb = "breeding_site_type", cloud = TRUE, ellipse = TRUE))
})
