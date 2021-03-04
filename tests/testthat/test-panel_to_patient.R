test_that("panel_to_patient works", {

  p2p <- panel_to_patient(gr_toy, mut_toy)

  expect_true(nrow(p2p) > 0)
  expect_true(ncol(p2p) > 0)
  expect_true(any(is.na(p2p$reg_id)))
})


test_that("panel_to_patient works with sv", {

  p2p <- panel_to_patient(gr_toy, mut_toy, sv_toy)


  expect_true(nrow(p2p) > 0)
  expect_true(ncol(p2p) > 0)
  expect_true("bp1_chr" %in% names(p2p))
  expect_true(any(is.na(p2p$reg_id)))
})

test_that("panel_to_patient works with sv and mode both", {

  p2p <- panel_to_patient(gr_toy, mut_toy, sv_toy)
  p2p_both <- panel_to_patient(gr_toy, mut_toy, sv_toy, sv_mode = "both")

  p2p_cov <- p2p %>%
    filter(!is.na(mut_id)) %>%
    pull(patient_id) %>%
    n_distinct()
  p2p_both_cov <- p2p_both %>%
    filter(!is.na(mut_id)) %>%
    pull(patient_id) %>%
    n_distinct()

  expect_true(p2p_cov >= p2p_both_cov)

})


