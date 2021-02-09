
#' Evaluate and annotate a given region to patient association data set
#'
#' @param reg_to_patient_mut A data.frame like object from `panel_to_patient()`
#'  with all regions `reg_id` linked to patients `patient_id` via mutations `mut_id`.
#'  It should additionally contain the following columns:
#'    - `reg_chr`
#'    - `reg_start` one-based
#'    - `reg_end`  one-based
#'
#' @return A data frame with all unique regions and additional annotations.
#'
#' @examples
#'
#' r2p <- panel_to_patient(gr_toy, mut_toy)
#' r2p
#'
#' eval_regions(r2p)
#'
#' @import dplyr
#' @export
eval_regions <- function(reg_to_patient_mut){

  # basic region annotation
  unique_reg <- reg_to_patient_mut %>%
    filter(!is.na(reg_id)) %>%
    distinct(reg_id, reg_chr, reg_start, reg_end) %>%
    mutate(
      n_region = row_number(),
      reg_size = reg_end - reg_start + 1,
      reg_size_cum = cumsum(reg_size)
    )

  # total number of patients for each region
  reg_to_n_patient <- reg_to_patient_mut %>%
    distinct(reg_id, patient_id) %>%
    count(reg_id, name = "n_patient")

  # total number of mutations for each region
  reg_to_n_mut <- reg_to_patient_mut %>%
    filter(!is.na(reg_id)) %>%
    distinct(reg_id, mut_id) %>%
    count(reg_id, name = "n_mut")

  # cumulative number of patietns for each region
  reg_to_n_patient_cum <- reg_to_patient_mut %>%
    filter(!is.na(reg_id)) %>%
    group_by(patient_id) %>%
    mutate(patient_counter = ifelse(row_number()==1,1,0)) %>%
    ungroup() %>%
    mutate(n_patients_cum = cumsum(patient_counter)) %>%
    group_by(reg_id) %>%
    summarize(
      n_patients_cum = max(n_patients_cum)
    )

  # cumulative number of mutations for each region
  reg_to_n_mut_cum <- reg_to_patient_mut %>%
    filter(!is.na(reg_id)) %>%
    group_by(mut_id) %>%
    mutate(mut_indicator = ifelse(row_number()==1,1,0)) %>%
    ungroup() %>%
    mutate(n_mut_cum = cumsum(mut_indicator)) %>%
    group_by(reg_id) %>%
    summarize(
      n_mut_cum = max(n_mut_cum)
    )

  reg_eval <- unique_reg %>%
    left_join(reg_to_n_patient, by = "reg_id") %>%
    left_join(reg_to_n_patient_cum, by = "reg_id") %>%
    mutate(
      total_patients = reg_to_patient_mut %>%
        pull(patient_id) %>%
        n_distinct(na.rm = TRUE),
      percent_patients = n_patient / total_patients * 100,
      percent_patients_cum = n_patients_cum / total_patients * 100,
      n_patient_gain = n_patients_cum - lag(n_patients_cum, default = 0),
    ) %>%
    left_join(reg_to_n_mut, by = "reg_id") %>%
    left_join(reg_to_n_mut_cum, by = "reg_id")

  return(reg_eval)
}

