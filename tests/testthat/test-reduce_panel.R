test_that("reduce_panel works with custom example data", {


  # example data ===============================================================
  #
  #    1 2 3 4 5 6 7 8 9 1011
  #r   |------|    |--|    ||
  #t     |------|    |--|
  #=============================================================================
  reg_gr <- GenomicRanges::GRanges(c(
    "1:100-500",
    "1:700-900",
    "1:1100-1200"
  ))
  target_gr <- GenomicRanges::GRanges(c(
    "1:200-600",
    "1:800-1000"
  ))

  r <- reduce_panel(reg_gr, target_gr, min_size = 100)
  expect_true(all(IRanges::overlapsAny(r, target_gr)))
  expect_true(all(S4Vectors::width(r) >= 100))

  r <- reduce_panel(reg_gr, target_gr, min_size = 200)
  expect_true(all(IRanges::overlapsAny(r, target_gr)))
  expect_true(all(S4Vectors::width(r) >= 200))

})
