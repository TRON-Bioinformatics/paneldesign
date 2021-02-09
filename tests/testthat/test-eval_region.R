test_that("eval_region works", {

  r2p <- panel_to_patient(gr_toy, mut_toy)
  r2p

  reg_eval <- eval_regions(r2p)

  expect_true(nrow(reg_eval) > 0)
  expect_true(ncol(reg_eval) > 0)
})
