
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
#' @param sv A data.frame like object with all structural variations and patients of the
#'   cohort of interest. Should have the following columns:
#'    - `patient_id`
#'    - `mut_id`
#'    - `bp1_chr`
#'    - `bp1_pos`
#'    - `bp2_chr`
#'    - `bp1_chr`
#'
#' @param sv_mode A character specifing if both brekapoints need to be coverd
#'  ("both") or a single breakpoint is sufficient ("single")
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
panel_to_patient <- function(reg_gr, mut, sv = NULL, sv_mode = "single"){

  # check input arguments
  stopifnot(sv_mode %in% c("both", "single"))

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

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Processing mutations
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Processing structural variants
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if(!is.null(sv)){
    # GRanges objec of all mutations
    bp1_gr <- GenomicRanges::GRanges(sv$bp1_chr,
                                     IRanges::IRanges(sv$bp1_pos, sv$bp1_pos),
                                     name = sv$mut_id)
    bp2_gr <- GenomicRanges::GRanges(sv$bp2_chr,
                                     IRanges::IRanges(sv$bp2_pos, sv$bp2_pos),
                                     name = sv$mut_id)

    # overlap regions with breakpoints
    reg_to_bp1 <- GenomicRanges::findOverlaps(reg_gr, bp1_gr) %>%
      as.data.frame() %>%
      dplyr::as_tibble() %>%
      dplyr::mutate(
        reg_id = names(reg_gr)[queryHits],
        mut_id = sv$mut_id[subjectHits],
      ) %>%
      dplyr::select(reg_id, mut_id)

    reg_to_bp2 <- GenomicRanges::findOverlaps(reg_gr, bp2_gr) %>%
      as.data.frame() %>%
      dplyr::as_tibble() %>%
      dplyr::mutate(
        reg_id = names(reg_gr)[queryHits],
        mut_id = bp2_gr$name[subjectHits],
      ) %>%
      dplyr::select(reg_id, mut_id)

    # combine breakpoint overlap according to sv mode
    if(sv_mode == "single"){
      # take union of region to mutation id associations
      reg_to_sv <- union(reg_to_bp1, reg_to_bp2)
    }else{
      # take intersection of region to mutation id associations
      reg_to_sv <- intersect(reg_to_bp1, reg_to_bp2)
    }

    # add regions not mapping to any sv
    reg_to_sv <- reg_df %>%
      dplyr::left_join(reg_to_sv, by = "reg_id")
  }else{
    reg_to_sv = NULL
  }
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Associate regions to patients
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # combine mutation and sv data
  mut_sv <- bind_rows(mut, sv)

  # mapping of region to pateint (that has at least one mutation in region)-----
  reg_to_patient <- dplyr::bind_rows(
      reg_to_mut,
      reg_to_sv
    ) %>%
    dplyr::left_join(
        mut_sv,
        by = "mut_id"
      )

  # add patients not mapping to any region
  reg_to_patient <- reg_to_patient %>%
    dplyr::full_join(
      dplyr::distinct(mut_sv, patient_id),
      by = "patient_id"
    )


  return(reg_to_patient)

}
