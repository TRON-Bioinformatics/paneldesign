

#' Create unique regions of fixed size around mutations as GRanges
#'
#' @param mut A data.frame like object with all mutations and patients of the
#'   cohort of interest. Should have the following columns:
#'    - `patient_id`
#'    - `mut_id`
#'    - `chr`
#'    - `start` or `pos`
#'    - `end` or `pos`
#' @param reg_size integer size of the region around mutations
#'
#' @examples
#' data(mut_toy)
#' mut_to_reg_gr(mut_toy)
#' mut_to_reg_gr(mut_toy, reg_size = 50)
#'
#' @import dplyr
#' @export
mut_to_reg_gr <- function(mut, reg_size = 100){


  stopifnot("start" %in% names(mut) | "pos" %in% names(mut) )
  stopifnot("end" %in% names(mut) | "pos" %in% names(mut) )

  # if pos is not given, use start
  if (! "pos" %in% names(mut)){
    mut$pos <- mut$start
  }

  # take unique mutation positions (by keeping Gene annotations) -----------------
  mut_pos <- mut %>%
    distinct(chr, pos) %>%
    mutate(
      mut_pos_id = row_number(),
      mut_pos_name = stringr::str_c(chr, pos, sep = "_")
    ) %>%
    select(mut_pos_id, chr, pos, mut_pos_name)

  # build Genomic ranges around unique mutations
  reg_gr <- GenomicRanges::GRanges(mut_pos$chr,
                                   IRanges::IRanges(
                                     mut_pos$pos - (reg_size/2 - 1),
                                     mut_pos$pos + reg_size/2
                                   ),
                                   name = mut_pos$mut_pos_name)

  return(reg_gr)
}

#' Create GRanges object out of mutations
#'
#' @param mut A data.frame like object with all mutations and patients of the
#'   cohort of interest. Should have the following columns:
#'    - `patient_id`
#'    - `mut_id`
#'    - `chr`
#'    - `start` or `pos`
#'    - `end` or `pos`
#' @examples
#' data(mut_toy)
#' mut_to_gr(mut_toy)
#'
#' @export
mut_to_gr <- function(mut){

  stopifnot("start" %in% names(mut) | "pos" %in% names(mut) )
  stopifnot("end" %in% names(mut) | "pos" %in% names(mut) )
  stopifnot("chr" %in% names(mut))
  stopifnot("mut_id" %in% names(mut))

  # if pos is not given, use start
  if (! "pos" %in% names(mut)){
    mut$pos <- mut$start
  }

  # GRanges objec of all mutations
  mut_gr <- GenomicRanges::GRanges(mut$chr,
                                   IRanges::IRanges(mut$pos, mut$pos),
                                   name = mut$mut_id)
  return(mut_gr)
}

