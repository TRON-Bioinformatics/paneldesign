#' Run mutation region definition and annotation
#'
#' @param mut a data.frame like object with mutations
#'
#' @import dplyr
#'
define_regions <- function(mut, reg_size = 100){

  # take unique mutation positions (by keeping Gene annotations) -----------------
  mut_pos <- mut %>%
    distinct(chr, pos, gene) %>%
    arrange(chr, pos, gene) %>%
    mutate(
      mut_pos_id = row_number(),
      mut_pos_name = stringr::str_c(gene, chr, pos, sep = "_")
    ) %>%
    select(mut_pos_id, chr, pos, gene, mut_pos_name)

  # build Genomic ranges around unique mutations
  reg_gr <- GenomicRanges::GRanges(mut_pos$chr,
                    IRanges::IRanges(
                      mut_pos$pos - (reg_size/2 - 1),
                      mut_pos$pos + reg_size/2
                    ),
                    name = mut_pos$mut_pos_name)

  # GRanges objec of all mutations
  mut_gr <- GenomicRanges::GRanges(mut$chr,
                    IRanges::IRanges(mut$pos, mut$pos),
                    name = mut$mut_id)

  # unique mutation region to all overlapping mutations --------------------------
  reg_to_mut <- GenomicRanges::findOverlaps(reg_gr, mut_gr) %>%
    as.data.frame() %>%
    as_tibble() %>%
    mutate(
      mut_pos_name = mut_pos$mut_pos_name[queryHits],
      mut_id = mut$mut_id[subjectHits],
    ) %>%
    select(mut_pos_name, mut_id)

  # mapping of region to pateint (that has at least one mutation in region)-----
  reg_to_patient <- reg_to_mut %>%
    left_join(mut, by = "mut_id") %>%
    distinct(cohort, mut_pos_name, patient_id)


  # Annotate regions ==========================================================


  # add number of patients with mutation in region
  reg_patient_count <- reg_to_patient %>%
    count(mut_pos_name, name = "n_reg_all_patient")

  # add number of mutations falling within this region
  reg_mut_count <- reg_to_mut %>%
    count(mut_pos_name, name = "n_reg_all_mut")

  # get region coordinates as character vector
  reg_to_coord <- tibble(
    mut_pos_name = reg_gr$name,
    reg_coord = as.character(reg_gr)
  )

  # get number of patients with mutated gene
  gene_to_patient_count <- mut_pos %>%
    left_join(reg_to_patient, by = "mut_pos_name") %>%
    distinct(gene, patient_id) %>%
    # filter Unknown genes out
    filter(gene != "Unknown") %>%
    count(gene, name = "patients_per_gene")

  # region to gene and patients with gene
  reg_to_gene <- mut_pos %>%
    distinct(mut_pos_name, gene) %>%
    left_join(gene_to_patient_count, by = "gene") %>%
    distinct(mut_pos_name, patients_per_gene)

  # total number of patient in all cohort
  n_all_patients <- mut %>%
    distinct(patient_id) %>%
    nrow()

  reg_annot <- mut_pos %>%
    left_join(reg_to_coord, by = "mut_pos_name") %>%
    left_join(reg_patient_count, by = "mut_pos_name") %>%
    left_join(reg_mut_count, by = "mut_pos_name") %>%
    left_join(reg_to_gene, by = "mut_pos_name") %>%
    # add total number of patients in all cohorts
    mutate(
      n_all_patients = n_all_patients
    )

  return(
    list(
      "reg_annot" = reg_annot,
      "reg_to_patient" = reg_to_patient
    )
  )

}
