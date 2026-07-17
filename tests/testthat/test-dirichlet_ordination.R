data(adults)
data(larvae)

test_that("dirichlet_ordination (adults) — title only across options", {
  withr::local_options(device = function(...) grDevices::png(filename = tempfile(fileext = ".png")))
  set.seed(20)

  EXPECTED_TITLE <- "Aligned Ordinations with Consensus Overlaid"

  # 1) group=location, color=location, cloud+ellipse
  p1 <- dirichlet_ordination(
    adults, draws = 10,
    group = "location", colorb = "location",
    cloud = TRUE, ellipse = TRUE
  )
  expect_true(!is.null(p1$plot), "result$plot is missing")
  expect_s3_class(p1$plot, "ggplot")
  expect_identical(p1$plot$labels$title, EXPECTED_TITLE)

  # 2) group=location, color=location, cloud only
  p2 <- dirichlet_ordination(
    adults, draws = 10,
    group = "location", colorb = "location",
    cloud = TRUE, ellipse = FALSE
  )
  expect_s3_class(p2$plot, "ggplot")
  expect_identical(p2$plot$labels$title, EXPECTED_TITLE)

  # 3) group=location, color=location, no cloud/ellipse
  p3 <- dirichlet_ordination(
    adults, draws = 10,
    group = "location", colorb = "location",
    cloud = FALSE, ellipse = FALSE
  )
  expect_s3_class(p3$plot, "ggplot")
  expect_identical(p3$plot$labels$title, EXPECTED_TITLE)

  # 4) group=location, color=sample_id, cloud only
  p4 <- dirichlet_ordination(
    adults, draws = 10,
    group = "location", colorb = "sample_id",
    cloud = TRUE, ellipse = FALSE
  )
  expect_s3_class(p4$plot, "ggplot")
  expect_identical(p4$plot$labels$title, EXPECTED_TITLE)

  # 5) group=location, color=sample_id, no cloud/ellipse
  p5 <- dirichlet_ordination(
    adults, draws = 10,
    group = "location", colorb = "sample_id",
    cloud = FALSE, ellipse = FALSE
  )
  expect_s3_class(p5$plot, "ggplot")
  expect_identical(p5$plot$labels$title, EXPECTED_TITLE)

  # 6) returned list has the expected structure
  expect_identical(p1$draws, 10)
  expect_true(is.data.frame(p1$df_consensus_coordinates))
  expect_true(is.data.frame(p1$df_all))
})

test_that("dirichlet_ordination (larvae) — different grouping column", {
  withr::local_options(device = function(...) grDevices::png(filename = tempfile(fileext = ".png")))
  set.seed(20)

  EXPECTED_TITLE <- "Aligned Ordinations with Consensus Overlaid"

  p6 <- dirichlet_ordination(
    larvae, draws = 10,
    group = "breeding_site_type", colorb = "breeding_site_type",
    cloud = TRUE, ellipse = TRUE
  )
  expect_true(!is.null(p6$plot), "result$plot is missing")
  expect_s3_class(p6$plot, "ggplot")
  expect_identical(p6$plot$labels$title, EXPECTED_TITLE)
})

test_that("dirichlet_ordination — too few draws to draw ellipse warns and disables it", {
  withr::local_options(device = function(...) grDevices::png(filename = tempfile(fileext = ".png")))
  set.seed(20)

  expect_warning(
    p7 <- dirichlet_ordination(
      adults, draws = 3,
      group = "location", colorb = "location",
      cloud = TRUE, ellipse = TRUE
    ),
    regexp = "Too few MC draws"
  )
  expect_s3_class(p7$plot, "ggplot")
})
