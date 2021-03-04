## code to prepare `sv_toy` dataset goes here


sv_toy <- tibble::tribble(
  ~patient_id, ~mut_id, ~bp1_chr, ~bp1_pos, ~bp2_chr, ~bp2_pos,
  "p4",        "sv01",  "1",      1000,    "1",      1500,
  "p5",        "sv02",  "1",      2000,    "2",      2000,
  )

usethis::use_data(sv_toy, overwrite = TRUE)
