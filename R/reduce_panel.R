
#' Reduce panel regions to input target ranges by respecting a minimum region size.
#'
#' @param reg_gr A GRanges object with regions defining the panel.
#'
#' @param target_gr A GRanges object with regions defining the target regions to
#' which the panel should be reduced.
#'
#' @param min_size integer size in bp as the minimal size an individual region
#'  can be reduced to.
#'
#' In cases where the overlapping part of an input region with a target region
#' is smaller than `min_size`, the overlapping part will be resized to `min_size`
#' by fixing the center of the overlapping part.
#'
#'
#' @examples
#' reg_gr <- GenomicRanges::GRanges(c(
#'   "1:100-500",
#'   "1:700-900",
#'   "1:1100-1200"
#' ))
#' target_gr <- GenomicRanges::GRanges(c(
#'   "1:200-600",
#'   "1:800-1000"
#' ))
#'
#' r <- reduce_panel(reg_gr, target_gr, min_size = 100)
#'
#' @export
reduce_panel <- function(reg_gr, target_gr, min_size = 100){

  # reduce target regiosn to non-overlapping regions
  target_gr_red <- GenomicRanges::reduce(target_gr, ignore.strand = TRUE)

  # compute pairs by any overlap
  pairs <- IRanges::findOverlapPairs(reg_gr, target_gr_red, ignore.strand = TRUE)

  # restrict to pairwise intersecting bases
  reduced <- IRanges::pintersect(pairs, ignore.strand = TRUE)

  # compute size of intersected regions
  size <- S4Vectors::width(reduced)

  # build results by
  #  - large enough intersected regions and
  #  - if smaller than min_size: resize overlapping regions to required min_size
  result <- c(
    reduced[size >= min_size],
    GenomicRanges::resize(reduced[size < min_size], min_size, fix = "center")
  )

  return(result)
}
