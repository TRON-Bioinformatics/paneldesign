test_that("select_greedy works", {

  mut_toy

  selection <- select_greedy(mut_toy)

  expect_equal(nrow(selection), 2)
  expect_equal(names(selection)[1], "mut_id")
  expect_true(all(selection$mut_id %in% mut_toy$mut_id))


})
