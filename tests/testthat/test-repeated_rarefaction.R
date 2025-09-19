data(adults)
data(larvae)

test_that("repeated_rarefaction (adults) — title only across options", {
  withr::local_options(device = function(...) grDevices::png(filename = tempfile(fileext = ".png")))
  set.seed(20)
  
  EXPECTED_TITLE <- "Aligned Ordinations with Consensus Overlaid"
  
  # 1) group=location, color=location, cloud+ellipse (RR_location_CLOUD_ELLIPSE)
  p1 <- repeated_rarefaction(
    adults, repeats = 10, threshold = 100,
    group = "location", colorb = "location",
    cloud = TRUE, ellipse = TRUE, seed = 20
  )
  expect_true(!is.null(p1$plot), "result$plot is missing")
  expect_s3_class(p1$plot, "ggplot")
  expect_identical(p1$plot$labels$title, EXPECTED_TITLE)
  
  # 2) group=location, color=location, cloud only (RR_location_CLOUD)
  p2 <- repeated_rarefaction(
    adults, repeats = 10, threshold = 100,
    group = "location", colorb = "location",
    cloud = TRUE, ellipse = FALSE, seed = 20
  )
  expect_s3_class(p2$plot, "ggplot")
  expect_identical(p2$plot$labels$title, EXPECTED_TITLE)
  
  # 3) group=location, color=location, no cloud/ellipse (RR_location)
  p3 <- repeated_rarefaction(
    adults, repeats = 10, threshold = 100,
    group = "location", colorb = "location",
    cloud = FALSE, ellipse = FALSE, seed = 20
  )
  expect_s3_class(p3$plot, "ggplot")
  expect_identical(p3$plot$labels$title, EXPECTED_TITLE)
  
  # 4) group=sample_id, color=location, cloud+ellipse (RR_CLOUD_ELLIPSE)
  p4 <- repeated_rarefaction(
    adults, repeats = 10, threshold = 100,
    group = "sample_id", colorb = "location",
    cloud = TRUE, ellipse = TRUE, seed = 20
  )
  expect_s3_class(p4$plot, "ggplot")
  expect_identical(p4$plot$labels$title, EXPECTED_TITLE)
  
  # 5) group=sample_id, color=location, cloud only (RR_sample/location_CLOUD)
  p5 <- repeated_rarefaction(
    adults, repeats = 10, threshold = 100,
    group = "sample_id", colorb = "location",
    cloud = TRUE, ellipse = FALSE, seed = 20
  )
  expect_s3_class(p5$plot, "ggplot")
  expect_identical(p5$plot$labels$title, EXPECTED_TITLE)
  
  # 6) group=sample_id, color=location, no cloud/ellipse (RR_sample/location)
  p6 <- repeated_rarefaction(
    adults, repeats = 10, threshold = 100,
    group = "sample_id", colorb = "location",
    cloud = FALSE, ellipse = FALSE, seed = 20
  )
  expect_s3_class(p6$plot, "ggplot")
  expect_identical(p6$plot$labels$title, EXPECTED_TITLE)
  
  # 7) group=location, color=sample_id, cloud only (RR_sample/sample_CLOUD)
  p7 <- repeated_rarefaction(
    adults, repeats = 10, threshold = 100,
    group = "location", colorb = "sample_id",
    cloud = TRUE, ellipse = FALSE, seed = 20
  )
  expect_s3_class(p7$plot, "ggplot")
  expect_identical(p7$plot$labels$title, EXPECTED_TITLE)
  
  # 8) group=location, color=sample_id, no cloud/ellipse (RR_sample/sample)
  p8 <- repeated_rarefaction(
    adults, repeats = 10, threshold = 100,
    group = "location", colorb = "sample_id",
    cloud = FALSE, ellipse = FALSE, seed = 20
  )
  expect_s3_class(p8$plot, "ggplot")
  expect_identical(p8$plot$labels$title, EXPECTED_TITLE)
})

# Test behavior when rarefaction threshold is not met on the larval dataset

test_that("repeated_rarefaction (larvae, threshold unmet) — warns and has correct title", {
  withr::local_options(device = function(...) grDevices::png(filename = tempfile(fileext = ".png")))
  set.seed(20)
  
  EXPECTED_TITLE <- "Aligned Ordinations with Consensus Overlaid"
  
  # RR_threshold_unmet
  expect_warning({
    p9 <- repeated_rarefaction(
      larvae, repeats = 10, threshold = 700,
      group = "breeding_site_type", colorb = "breeding_site_type",
      cloud = TRUE, ellipse = TRUE, seed = 20
    )
    expect_true(!is.null(p9$plot), "result$plot is missing")
    expect_s3_class(p9$plot, "ggplot")
    expect_identical(p9$plot$labels$title, EXPECTED_TITLE)
  })
})
