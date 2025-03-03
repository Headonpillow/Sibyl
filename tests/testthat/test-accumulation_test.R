data(adults)
test_that("Acculumation_works", {
  vdiffr::expect_doppelganger("Accumulation plot", accumulation_test(adults, step = 50))
})

