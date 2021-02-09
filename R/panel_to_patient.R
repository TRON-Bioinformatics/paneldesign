
#' Get mapping of panel regions to covered patients and mutations.
#'
#' @param reg_gr A GRanges object with regions defining the panel. Should have
#'   unique names and is expected to be sorted by decreasing priority.
#' @param mut A data.frame like object with all mutations and patients of the
#'   cohort of interest. Should have the following columns:
#'    - `patient_id`
#'    - `mut_id`
#'    - `chr`
#'    - `start`
#'    - `end`
#'
#' @return A data set mapping the regions to affected patients and mutations.
#'  The output contains at least the following columns
#'
#'  - `reg_id`
#'  - `mut_id`
#'  - `patient_id`
#'
#' @examples
#' gr_toy
#' mut_toy
#'
#' panel_to_patient(gr_toy, mut_toy)
#'
#' @export
panel_to_patient <- function(reg_gr, mut){

  # if no names exist in region use running numbers
  if (is.null(names(reg_gr))){
    names(reg_gr) <- stringr::str_c("reg_", seq_along(reg_gr))
  }

  # assume names of regions to be unique
  stopifnot(
    length(unique(names(reg_gr))) == length(reg_gr)
  )

  # if start and end is not given, use pos
  if (! "start" %in% names(mut) & ! "end" %in% names(mut)){
    mut$start <- mut$pos
    mut$end <- mut$pos
  }

  reg_df <- tibble::tibble(
    reg_id = names(reg_gr),
    reg_chr = GenomeInfoDb::seqnames(reg_gr) %>% as.character(),
    reg_start = start(reg_gr),
    reg_end = end(reg_gr)
  )

  # GRanges objec of all mutations
  mut_gr <- GenomicRanges::GRanges(mut$chr,
                                   IRanges::IRanges(mut$start, mut$end),
                                   name = mut$mut_id)

  # overlap regions with mutations
  reg_to_mut <- GenomicRanges::findOverlaps(reg_gr, mut_gr) %>%
    as.data.frame() %>%
    dplyr::as_tibble() %>%
    dplyr::mutate(
      reg_id = names(reg_gr)[queryHits],
      mut_id = mut$mut_id[subjectHits],
    ) %>%
    dplyr::select(reg_id, mut_id)

  # add regions not mapping to any mutation
  reg_to_mut <- reg_df %>%
    dplyr::left_join(reg_to_mut, by = "reg_id")

  # mapping of region to pateint (that has at least one mutation in region)-----
  reg_to_patient <- reg_to_mut %>%
    dplyr::left_join(mut, by = "mut_id")

  # add patients not mapping to any region
  reg_to_patient <- reg_to_patient %>%
    dplyr::full_join(
      dplyr::distinct(mut, patient_id),
      by = "patient_id"
    )

  return(reg_to_patient)

}
