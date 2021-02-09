
#' Evaluate and annotate a given set of regions (panel)
#'
#' @param reg_to_patient_mut A data.frame like object from `panel_to_patient()`
#'  with all regions `reg_id` linked to patients `patient_id` via mutations `mut_id`.
#'  It should additionally contain the following columns:
#'    - `reg_chr`
#'    - `reg_start` one-based
#'    - `reg_end`  one-based
#'
#' @return A data frame with a single row and annotation columns to the panel
#'
#' @examples
#'
#' r2p <- panel_to_patient(gr_toy, mut_toy)
#' r2p
#'
#' eval_panel(r2p)
#'
#' @import dplyr
#' @export
eval_panel <- function(reg_to_patient_mut){

  mut_per_patient_df <- reg_to_patient_mut %>%
    filter(!is.na(patient_id)) %>%
    distinct(patient_id, mut_id) %>%
    group_by(patient_id) %>%
    summarize(
      mut_per_patient_n = sum(!is.na(mut_id))
    )

  tibble::tibble(
    n_reg = reg_to_patient_mut %>% pull(reg_id) %>% n_distinct(na.rm = TRUE),

    size_total = reg_to_patient_mut %>%
      filter(!is.na(reg_id)) %>%
      distinct(reg_id, reg_end, reg_start) %>%
      mutate(reg_size = reg_end - reg_start + 1) %>%
      pull(reg_size) %>% sum(),

    total_patients = reg_to_patient_mut %>%
      pull(patient_id) %>%
      n_distinct(na.rm = TRUE),

    n_patient = reg_to_patient_mut %>%
      filter(!is.na(reg_id)) %>%
      pull(patient_id) %>%
      n_distinct(na.rm = TRUE),

    percent_patient = n_patient / total_patients * 100,

    n_mut = reg_to_patient_mut %>% pull(mut_id) %>% n_distinct(na.rm = TRUE),

    mut_per_patient_median = mut_per_patient_df %>%
      pull(mut_per_patient_n) %>%
      median(),

    mut_per_patient_mean = mut_per_patient_df %>%
      pull(mut_per_patient_n) %>%
      mean(),

    mut_per_patient_df = list(mut_per_patient_df),
  )


}
