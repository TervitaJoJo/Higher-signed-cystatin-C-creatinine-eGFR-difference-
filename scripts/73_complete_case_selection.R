# Missingness-selection inverse probability weighting and HRS version/overlap audit.

suppressPackageStartupMessages({
  library(dplyr)
  library(survey)
  library(jsonlite)
})
options(survey.lonely.psu = "adjust", survey.adjust.domain.lonely = TRUE)

root <- normalizePath(commandArgs(TRUE)[1], winslash="/", mustWork=TRUE)
dir.create(file.path(root, "audit"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(root, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(root, "reports"), recursive = TRUE, showWarnings = FALSE)

cache <- readRDS(file.path(root, "data_internal", "common_endpoint_models.rds"))

make_design <- function(d, cohort) {
  if (startsWith(cohort, "HRS")) {
    svydesign(ids = ~psu, strata = ~stratum, weights = ~ipw_weight,
              data = d, nest = TRUE)
  } else {
    svydesign(ids = ~psu, weights = ~ipw_weight, data = d)
  }
}

fit_svy <- function(d, cohort) {
  predictors <- c("delta10", "mean_egfr", "age_c", "female", "education",
                  "interval", "start_severity", "bmi", "smoking",
                  "hypertension", "diabetes", "heart", "stroke")
  if (cohort == "HRS2016VBS") predictors <- c(predictors, "draw_lag_months")
  z <- d[complete.cases(d[, c("y", predictors, "psu", "ipw_weight",
                              if (startsWith(cohort, "HRS")) "stratum"), drop = FALSE]), , drop = FALSE]
  z <- z[is.finite(z$ipw_weight) & z$ipw_weight > 0, , drop = FALSE]
  z$ipw_weight <- z$ipw_weight / mean(z$ipw_weight)
  design <- make_design(z, cohort)
  fit <- tryCatch(
    svyglm(reformulate(predictors, response = "y"), design = design,
           family = quasibinomial()),
    error = function(e) e
  )
  if (inherits(fit, "error")) {
    return(data.frame(cohort = cohort, model = NA_character_, n_people = length(unique(z$pid)),
                      n_intervals = nrow(z), events = sum(z$y == 1), odds_ratio = NA,
                      ci_lower = NA, ci_upper = NA, p = NA, stringsAsFactors = FALSE))
  }
  cf <- summary(fit)$coefficients["delta10", ]
  ci <- suppressMessages(confint(fit, parm = "delta10"))
  data.frame(cohort = cohort, model = NA_character_, n_people = length(unique(z$pid)),
             n_intervals = nrow(z), events = sum(z$y == 1),
             odds_ratio = exp(unname(cf["Estimate"])), ci_lower = exp(ci[1]),
             ci_upper = exp(ci[2]), p = unname(cf["Pr(>|t|)"]), stringsAsFactors = FALSE)
}

prepare_selection <- function(d, vars) {
  s <- d
  for (v in vars) {
    if (is.factor(s[[v]]) || is.character(s[[v]])) {
      value <- as.character(s[[v]])
      value[is.na(value) | !nzchar(value)] <- "__MISS__"
      s[[v]] <- factor(value)
    } else {
      miss <- is.na(s[[v]]) | !is.finite(s[[v]])
      s[[paste0(v, "_miss")]] <- as.integer(miss)
      if (any(miss)) s[[v]][miss] <- median(s[[v]][!miss], na.rm = TRUE)
    }
  }
  s
}

run_ipw <- function(key, cohort) {
  d <- cache[[key]][["data"]]
  d <- d[d$start_state == "difficulty", , drop = FALSE]
  d$delta <- d$egfr_cys - d$egfr_cr
  d$delta10 <- d$delta / 10
  d$mean_egfr <- (d$egfr_cys + d$egfr_cr) / 2
  d$y <- as.integer(d$next_state == "independent")
  d$known_outcome <- !is.na(d$next_state) & d$next_state != "unknown"

  predictors <- c("delta10", "mean_egfr", "age_c", "female", "education",
                  "interval", "start_severity", "bmi", "smoking",
                  "hypertension", "diabetes", "heart", "stroke")
  if (cohort == "HRS2016VBS") predictors <- c(predictors, "draw_lag_months")
  required <- c("y", predictors, "psu", "analysis_weight")
  if (startsWith(cohort, "HRS")) required <- c(required, "stratum")
  d$complete_analysis <- complete.cases(d[, required, drop = FALSE]) & d$known_outcome

  sel_vars <- c("age_c", "female", "education", "bmi", "smoking",
                "hypertension", "diabetes", "heart", "stroke",
                "start_severity", "interval")
  if (cohort == "HRS2016VBS") sel_vars <- c(sel_vars, "draw_lag_months")
  eligible <- d$known_outcome & is.finite(d$analysis_weight) & d$analysis_weight > 0
  s <- prepare_selection(d[eligible, , drop = FALSE], sel_vars)
  sel_terms <- c(sel_vars, paste0(sel_vars[vapply(sel_vars, function(v) is.numeric(d[[v]]), logical(1))], "_miss"))
  sel_terms <- unique(sel_terms[sel_terms %in% names(s)])
  sel_formula <- reformulate(sel_terms, response = "complete_analysis")
  sel_fit <- tryCatch(glm(sel_formula, data = s, family = binomial()), error = function(e) e)
  if (inherits(sel_fit, "error")) stop("Selection model failed for ", cohort, ": ", conditionMessage(sel_fit))

  s$p_complete <- as.numeric(predict(sel_fit, newdata = s, type = "response"))
  s$p_complete <- pmin(pmax(s$p_complete, 0.05), 0.95)
  d$p_complete <- NA_real_
  d$p_complete[eligible] <- s$p_complete
  z <- d[d$complete_analysis, , drop = FALSE]
  z$p_complete <- d$p_complete[d$complete_analysis]
  z$ipw_weight <- z$analysis_weight * mean(z$p_complete, na.rm = TRUE) / z$p_complete
  q <- quantile(z$ipw_weight, c(0.01, 0.99), na.rm = TRUE)
  z$ipw_weight <- pmin(pmax(z$ipw_weight, q[1]), q[2])
  ipw <- fit_svy(z, cohort)
  ipw$model <- "selection_IPW_trim1_99"
  z$ipw_weight <- z$analysis_weight
  primary <- fit_svy(z, cohort)
  primary$model <- "complete_case_reference"

  out <- rbind(primary, ipw)
  out$eligible_known_intervals <- sum(eligible)
  out$complete_intervals <- sum(d$complete_analysis)
  out$complete_fraction <- mean(d$complete_analysis[eligible])
  out$selection_p_median <- median(s$p_complete, na.rm = TRUE)
  out$selection_p_p01 <- unname(quantile(s$p_complete, 0.01, na.rm = TRUE))
  out$selection_p_p99 <- unname(quantile(s$p_complete, 0.99, na.rm = TRUE))
  out$ipw_q01 <- unname(q[1])
  out$ipw_q99 <- unname(q[2])
  out$ipw_effective_n <- (sum(z$ipw_weight)^2) / sum(z$ipw_weight^2)
  list(result = out, selection = s, data = z)
}

runs <- list(
  CHARLS2011 = run_ipw("CHARLS2011_egfr_cr", "CHARLS2011"),
  HRS2016VBS = run_ipw("HRS2016VBS_egfr_cr", "HRS2016VBS")
)
ipw_results <- bind_rows(lapply(runs, `[[`, "result"))
write.csv(ipw_results, file.path(root, "tables", "discordance_missingness_ipw_results.csv"), row.names = FALSE)

