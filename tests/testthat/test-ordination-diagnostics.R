test_that("ordination preserves native geometry and reports reference eigenvalues", {
  counts <- rbind(a = c(12, 3, 1, 4), b = c(2, 14, 3, 1),
                  c = c(1, 2, 15, 2), d = c(5, 1, 2, 12))
  out <- Sibyl:::ord_and_mean(list(counts, counts), 2, cores = 1)
  expected <- stats::cmdscale(vegan::vegdist(counts), k = 2, eig = TRUE)
  expect_equal(as.vector(stats::dist(out$consensus_coordinates)),
               as.vector(stats::dist(expected$points)))
  expect_equal(out$axis_variance,
               100 * pmax(expected$eig[1:2], 0) / sum(abs(expected$eig)))
  info <- data.frame(location = c("B", "A", "B", "A"),
                     row.names = c("d", "c", "b", "a"))
  plotted <- Sibyl:::plot_rep_raref(out$aligned_ordinations,
    out$consensus_coordinates, info, "location", "location",
    TRUE, FALSE, "test", out$axis_variance)
  expect_equal(plotted$consensus_df$location,
               info[plotted$consensus_df$sample_id, "location"])
  expect_match(plotted$plot$labels$x, "%", fixed = TRUE)
  expect_equal(plotted$plot$coordinates$ratio, 1)
})

test_that("CH uses native consensus geometry and excluded APD is missing", {
  withr::local_options(device = function(...) grDevices::png(filename = tempfile(fileext = ".png")))
  counts <- matrix(c(12,3,1,4, 2,14,3,1, 1,2,15,2, 5,1,2,12,
                     1,1,1,1), nrow = 5, byrow = TRUE,
                   dimnames = list(letters[1:5], paste0("tax", 1:4)))
  info <- data.frame(location = c("A", "A", "B", "B", "A"),
                     row.names = letters[1:5])
  ps <- phyloseq::phyloseq(phyloseq::otu_table(t(counts), taxa_are_rows = TRUE),
                          phyloseq::sample_data(info))
  set.seed(27)
  expect_warning(rr <- repeated_rarefaction(ps, repeats = 5, threshold = 10,
       group = "location", colorb = "location", cores = 1), "removed")
  expected <- clusterSim::index.G1(rr$df_consensus_coordinates[, c("Dim1", "Dim2")],
                                  c(1, 1, 2, 2))
  set.seed(27)
  expect_warning(scan <- test_threshold(ps, repeats = 5,
       t_min = 10, t_max = 10, group = "location", cores = 1, verbose = FALSE),
       "removed")
  expect_named(scan, c("index_plot", "ordination_plots", "avg_distances", "index_data"))
  expect_equal(scan$index_data$Index, expected)
  apd <- scan$avg_distances$repeat_number_5
  expect_true(is.na(apd$Avg_Distance[apd$Sample_ID == "e"]))
  expect_true(all(is.finite(apd$Avg_Distance[apd$Sample_ID != "e"])))
  expect_equal(avg_pairwise_dist_plot(apd)$theme$legend.position, "none")
  expect_equal(avg_pairwise_dist_plot(apd, show_legend = TRUE)$theme$legend.position, "right")
})

test_that("count extraction respects either phyloseq table orientation", {
  data("adults", package = "Sibyl", envir = environment())
  ps <- adults
  expected <- methods::as(phyloseq::otu_table(ps), "matrix")
  if (phyloseq::taxa_are_rows(ps)) expected <- t(expected)
  expect_equal(Sibyl:::sample_count_matrix(ps), expected)
  phyloseq::otu_table(ps) <- phyloseq::otu_table(expected, taxa_are_rows = FALSE)
  expect_equal(Sibyl:::sample_count_matrix(ps), expected)
})

test_that("seeded replicates and native ordinations agree across worker counts", {
  counts <- rbind(a = c(12,3,1,4), b = c(2,14,3,1),
                  c = c(1,2,15,2), d = c(5,1,2,12))
  set.seed(81)
  one <- Sibyl:::rep_raref(counts, threshold = 10, repeats = 3, cores = 1)
  set.seed(81)
  two <- Sibyl:::rep_raref(counts, threshold = 10, repeats = 3, cores = 2)
  expect_equal(one$rarefied_matrix_list, two$rarefied_matrix_list)
  a <- Sibyl:::ord_and_mean(one$rarefied_matrix_list, 3, cores = 1)
  b <- Sibyl:::ord_and_mean(two$rarefied_matrix_list, 3, cores = 2)
  expect_equal(a$consensus_coordinates, b$consensus_coordinates)
  expect_equal(a$axis_variance, b$axis_variance)
})
