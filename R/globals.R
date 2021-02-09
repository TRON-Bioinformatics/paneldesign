
# define global variables just to avoid notes in R CMD check
# See https://community.rstudio.com/t/how-to-solve-no-visible-binding-for-global-variable-note/28887

utils::globalVariables(c("end", "median", "mut_id", "mut_indicator",
                         "mut_per_patient_n", "n_cum", "n_mut_cum",
                         "n_patient", "n_patients_cum", "n_samples",
                         "patient_counter", "patient_id", "queryHits",
                         "reg_chr", "reg_end", "reg_id", "reg_size",
                         "reg_start", "start", "subjectHits", "total_patients"))

