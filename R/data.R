#' A toy example data set linking patients to mutations.
#'
#'
#' @format A data frame with rows and 5 columns:
#' \describe{
#'   \item{patient_id}{patient ID}
#'   \item{mut_id}{unique mutation ID}
#'   \item{chr}{chromosome name}
#'   \item{start}{start coordinate of the mutation}
#'   \item{end}{end coordinate of the mutation}
#' }
#' @source made up
"mut_toy"

#' A toy example genomic ranges object.
#'
#' This example *panel* contains regions overlapping toy example mutations `m1`
#' and `m2` from \code{\link{mut_toy}}.
#'
#'
#' @format A \code{\link[GenomicRanges]{GRanges}} object with rows 2 ranges.
#'
#' @source made up
"gr_toy"
