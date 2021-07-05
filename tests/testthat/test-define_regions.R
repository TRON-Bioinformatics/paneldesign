test_that("get_exons works", {

  # get hg19 ensembl annotations
  require(EnsDb.Hsapiens.v75)

  edb <- EnsDb.Hsapiens.v75

  # ex_full <- get_exons(edb)

  # query all exons from chr22
  ex_chr22 <- get_exons(
    edb,
    filter = AnnotationFilter::GRangesFilter(
      GenomicRanges::GRanges("22", IRanges(1, 51304566))
      )
    )

  expect_equal(as.character(class(ex_chr22)), "GRanges")

  # test that all are on chr22
  expect_equal(unique(GenomeInfoDb::seqlevels(ex_chr22)), "chr22")


})
