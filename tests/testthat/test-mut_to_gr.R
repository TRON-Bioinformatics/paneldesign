test_that("mut_to_gr works", {

  gr <- mut_to_gr(mut_toy)

  expect_equal(length(gr), nrow(mut_toy))
  expect_equal(IRanges::width(gr), rep(1, length(gr)))
})
