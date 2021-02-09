## code to prepare `mut_toy` dataset goes here


mut_toy <- tibble::tribble(
  ~patient_id, ~mut_id, ~chr, ~start, ~end,
  "p1",        "m01",   "1",  1000,   1000,
  "p1",        "m02",   "1",  2000,   2000,
  "p1",        "m03",   "2",  3000,   3000,

  "p2",        "m01",   "1",  1000,   1000,
  "p2",        "m04",   "2",  4000,   4000,
  "p2",        "m05",   "3",  5000,   5000,

  "p3",        "m02",   "1",  2000,   2000,
  "p3",        "m04",   "2",  4000,   4000,
  "p3",        "m06",   "2",  6000,   6000,
)

usethis::use_data(mut_toy, overwrite = TRUE)
