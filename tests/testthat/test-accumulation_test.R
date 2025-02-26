data(adults)
test_that("repeated_rarefaction_works", {
  vdiffr::expect_doppelganger("Accumulation plot", accumulation_test(adults, step = 50))
})

