test_that("panel_to_patient works", {

  p2p <- panel_to_patient(gr_toy, mut_toy)

  expect_true(nrow(p2p) > 0)
  expect_true(ncol(p2p) > 0)
  expect_true(any(is.na(p2p$reg_id)))
})
