test_that("mut_to_reg_gr works", {

  reg_gr <- mut_to_reg_gr(mut_toy, reg_size = 100)

  expect_lte(length(reg_gr), nrow(mut_toy))
  expect_equal(IRanges::width(reg_gr), rep(100, length(reg_gr)))

})
