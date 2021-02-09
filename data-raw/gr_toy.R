## code to prepare `gr_toy` dataset goes here

gr_toy <- GenomicRanges::GRanges(
  c(
    "1:1000-1000",
    "2:3000-3000"
    )
)
names(gr_toy) <- c(
  "r1",
  "r3"
)


usethis::use_data(gr_toy, overwrite = TRUE)
