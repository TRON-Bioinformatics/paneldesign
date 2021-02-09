

#' Select mutations by maximal coverage
#'
#' @param set_df A dataset of mutations (sets) in patients. Each line is an
#'           occurrence of a mutation in a patient.
#' @param mut_col Character for column name of the mutation identifier (set)
#' @param sample_col Character for column name of the sample identifier (elements in set)
#'
#' @return a [tibble][tibble::tibble-package] with the selected mutations (sets), number and coverage of the mutations
#'
#' @examples
#' mut_toy
#'
#' select_greedy(mut_toy)
#'
#' @export
select_greedy <- function(set_df, mut_col = "mut_id",
                          sample_col = "patient_id"){

  mut_set_cover <- set_df %>%
    dplyr::select(dplyr::all_of(c(mut_col, sample_col))) %>%
    RcppGreedySetCover::greedySetCover(FALSE) %>%
    dplyr::as_tibble()

  mut_set_cover %>%
    dplyr::group_by_at(mut_col) %>%
    dplyr::summarize(
      n = dplyr::n()
    ) %>%
    dplyr::arrange(dplyr::desc(n)) %>%
    dplyr::mutate(
      n_samples = sum(n),
      rank = rank(dplyr::desc(n)),
      order = order(n, decreasing = TRUE),
      n_cum = cumsum(n),
      coverage = n / n_samples * 100,
      coverage_cum = n_cum / n_samples * 100,
    )
}

