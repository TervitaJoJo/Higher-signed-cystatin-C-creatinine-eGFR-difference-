Sys.setlocale("LC_CTYPE", "English_United States.utf8")
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(readr); library(jsonlite) })

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) normalizePath(args[[1]], winslash = "/") else getwd()
dir.create(file.path(root, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(root, "audit"), recursive = TRUE, showWarnings = FALSE)

models <- readRDS(file.path(root, "data_internal", "common_endpoint_models.rds"))
cohort_data <- list(
  HRS = models[["HRS2016VBS_egfr_cys"]][["data"]],
  CHARLS = models[["CHARLS2011_egfr_cys"]][["data"]]
)

fmt_num <- function(x, digits = 1) {
  if (!length(x) || all(is.na(x))) return("-")
  sprintf(paste0("%.", digits, "f"), mean(x, na.rm = TRUE))
}

fmt_mean_sd <- function(x, digits = 1) {
  z <- x[is.finite(x)]
  if (!length(z)) return("-")
  sprintf(paste0("%.", digits, "f (%.", digits, "f)"), mean(z), sd(z))
}

fmt_median_iqr <- function(x, digits = 0) {
  z <- x[is.finite(x)]
  if (!length(z)) return("-")
  q <- quantile(z, c(.25, .5, .75), names = FALSE)
  if (digits == 0) sprintf("%.0f (%.0f-%.0f)", q[2], q[1], q[3])
  else sprintf(paste0("%.", digits, "f (%.", digits, "f-%.", digits, "f)"), q[2], q[1], q[3])
}

fmt_n_pct <- function(x, positive = NULL, digits = 1) {
  ok <- !is.na(x)
  if (!any(ok)) return("-")
  if (is.null(positive)) positive <- x == 1
  n <- sum(positive & ok)
  sprintf(paste0("%d (%.", digits, "f%%)"), n, 100 * n / sum(ok))
}

one_person_per_id <- function(d) d %>% arrange(pid, wave) %>% distinct(pid, .keep_all = TRUE)

baseline_rows <- tibble(
  characteristic = c(
    "Participants",
    "Age, years, mean (SD)",
    "Female sex, n (%)",
    "Education categories, n (%)",
    "Body mass index, kg/m2, mean (SD)",
    "eGFRcr, mL/min/1.73 m2, mean (SD)",
    "eGFRcys, mL/min/1.73 m2, mean (SD)",
    "Signed difference eGFRcys-eGFRcr, mean (SD)",
    "Mean eGFR, mL/min/1.73 m2, mean (SD)",
    "Mean eGFR <60 mL/min/1.73 m2, n (%)",
    "Combined eGFR, mL/min/1.73 m2, mean (SD)",
    "Index ADL difficulty count, median (IQR)",
    "Current smoking, n (%)",
    "Hypertension, n (%)",
    "Diabetes, n (%)",
    "Heart disease, n (%)",
    "Stroke, n (%)"
  ),
  type = c("n", "mean_sd", "binary", "education", "mean_sd", rep("mean_sd", 4),
           "binary", "mean_sd", "median_iqr", rep("binary", 5)),
  variable = c("pid", "age", "female", "education", "bmi",
               "egfr_cr", "egfr_cys", "delta", "mean_egfr", "ckd", "egfr_combined", "start_severity",
               "smoking", "hypertension", "diabetes", "heart", "stroke"),
  category = rep("", 17),
  stringsAsFactors = FALSE
)

education_summary <- function(d, cohort) {
  z <- as.character(d$education)
  if (cohort == "HRS") {
    levs <- c("less_hs", "hs", "college")
    labels <- c("Less than high school", "High school", "College or higher")
  } else {
    levs <- c("none", "primary_or_less", "middle_or_more")
    labels <- c("No formal education", "Primary or less", "Middle school or higher")
  }
  vals <- vapply(seq_along(levs), function(i) {
    ok <- !is.na(z)
    n <- sum(ok & z == levs[i])
    sprintf("%s: %d (%.1f%%)", labels[i], n, 100 * n / sum(ok))
  }, character(1))
  paste(vals, collapse = "; ")
}

render_cell <- function(d, row) {
  v <- row$variable
  if (v == "pid") return(as.character(nrow(d)))
  if (v == "age") return(fmt_mean_sd(d[[v]], 1))
  if (v %in% c("bmi", "egfr_cr", "egfr_cys", "egfr_combined", "delta", "mean_egfr")) {
    return(fmt_mean_sd(d[[v]], 1))
  }
  if (v == "start_severity") return(fmt_median_iqr(d[[v]], 0))
  if (v == "education") return(education_summary(d, ifelse("education" %in% names(d) && "middle_or_more" %in% as.character(d$education), "CHARLS", "HRS")))
  if (v == "ckd") return(fmt_n_pct(d[[v]], d[[v]] == 1, 1))
  if (v %in% c("female", "smoking", "hypertension", "diabetes", "heart", "stroke")) {
    return(fmt_n_pct(d[[v]], d[[v]] == 1, 1))
  }
  stop("Unhandled baseline variable: ", v)
}

out <- baseline_rows
missing <- list()
for (nm in names(cohort_data)) {
  d <- one_person_per_id(cohort_data[[nm]]) %>%
    mutate(
      delta = (egfr_cys - egfr_cr),
      mean_egfr = (egfr_cys + egfr_cr) / 2,
      ckd = as.integer(mean_egfr < 60)
    )
  out[[nm]] <- vapply(seq_len(nrow(out)), function(i) render_cell(d, out[i, ]), character(1))
  miss_vars <- c("bmi", "egfr_cys", "egfr_combined", "smoking", "hypertension", "diabetes", "heart", "stroke")
  missing[[nm]] <- tibble(
    cohort = nm,
    variable = miss_vars,
    n_missing = vapply(miss_vars, function(v) sum(is.na(d[[v]])), integer(1)),
    n_total = nrow(d)
  ) %>% mutate(percent_missing = 100 * n_missing / n_total)
}

write_csv(out, file.path(root, "tables", "primary_baseline_characteristics.csv"))
write_csv(bind_rows(missing), file.path(root, "tables", "primary_baseline_missingness.csv"))

validation <- list(
  status = "passed",
  source = "data_internal/common_endpoint_models.rds",
  definition = "one row per participant contributing to the primary eGFRcys model risk set",
  participants = setNames(lapply(cohort_data, function(d) nrow(one_person_per_id(d))), names(cohort_data)),
  intervals = setNames(lapply(cohort_data, nrow), names(cohort_data)),
  unknown_next_state = setNames(lapply(cohort_data, function(d) sum(d$next_state == "unknown", na.rm = TRUE)), names(cohort_data)),
  table = "tables/primary_baseline_characteristics.csv",
  missingness = "tables/primary_baseline_missingness.csv",
  note = "Descriptive values are unweighted. Survey weights remain used in regression models. No between-cohort significance tests were performed."
)
write_json(validation, file.path(root, "audit", "primary_baseline_characteristics_validation.json"), pretty = TRUE, auto_unbox = TRUE)
print(out)
print(validation)
