test_that("eval_panel works", {

  r2p <- panel_to_patient(gr_toy, mut_toy)
  r2p

  panel_eval <- eval_panel(r2p)

  expect_true(nrow(panel_eval) > 0)
  expect_true(ncol(panel_eval) > 0)
})
