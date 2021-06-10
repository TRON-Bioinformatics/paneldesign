## code to prepare `mut_toy` dataset goes here


mut_toy <- tibble::tribble(
  ~patient_id, ~mut_id, ~chr, ~start, ~end, ~gene,
  "p1",        "m01",   "1",  1000,   1000, "g1",
  "p1",        "m02",   "1",  2000,   2000, "g1",
  "p1",        "m03",   "2",  3000,   3000, "g2",

  "p2",        "m01",   "1",  1000,   1000, "g1",
  "p2",        "m04",   "2",  4000,   4000, "g2",
  "p2",        "m05",   "3",  5000,   5000, "g3",

  "p3",        "m02",   "1",  2000,   2000, "g1",
  "p3",        "m04",   "2",  4000,   4000, "g2",
  "p3",        "m06",   "2",  6000,   6000, "g2",
)

usethis::use_data(mut_toy, overwrite = TRUE)
