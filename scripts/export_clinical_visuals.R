Sys.setlocale("LC_CTYPE", "English_United States.utf8")
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(survey); library(splines); library(jsonlite)
})
options(survey.lonely.psu = "adjust", survey.adjust.domain.lonely = TRUE)
project <- normalizePath(Sys.getenv("REPRO_WORK"), winslash = "/")
root <- file.path(project, "analysis")
out <- file.path(root, "audit/bmc_figure_restructure_20261011/source_data")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

# Reuse the original cohort construction and model functions without running
# its result-writing loop or replacing any locked analysis output.
src <- readLines(file.path(root, "code/102_clinical_nonlinear_interaction_competing_risk.R"))
first <- grep('^cc <- readRDS', src)[1]
last <- grep('^model_rows <- list', src)[1] - 1L
eval(parse(text = src[first:last]))
locked <- read.csv(file.path(root, "tables/clinical_enhancement_model_results.csv"))
nonlinear <- read.csv(file.path(root, "tables/clinical_nonlinear_spline_results.csv"))
curves <- list(); samples <- list(); flows <- list(); quantiles <- list(); checks <- list()

for (cohort in c("HRS2016VBS", "CHARLS2011")) {
  co <- cohorts[[cohort]]
  raw <- intervals_for(co) %>% filter(start_state == "difficulty")
  d <- prepare(co, cohort)
  bt <- base_terms(d, cohort)
  f <- fit_svy(d, cohort, "recovery_event", bt, "primary_linear")
  z <- f$data
  target <- locked %>% filter(.data$cohort == .env$cohort, analysis == "primary_linear",
    outcome == "recovery_event", term == "delta10")
  stopifnot(nrow(z) == target$n_intervals, n_distinct(z$pid) == target$n_people,
    sum(z$recovery_event) == target$events,
    abs(coef(f$fit)["delta10"] - target$estimate) < 1e-10)

  # This is exactly the previously specified four-df spline and orthogonalisation.
  b <- ns(z$delta10, df = 4)
  residual_model <- lm(b ~ z$delta10)
  br <- residuals(residual_model)
  q <- qr(br)
  k <- min(3L, q$rank)
  nl <- qr.Q(q)[, seq_len(k), drop = FALSE]
  nl_names <- paste0("rcs_nl", seq_len(k))
  colnames(nl) <- nl_names
  z[, nl_names] <- nl
  sf <- fit_svy(z, cohort, "recovery_event", c(bt, nl_names), "natural_cubic_spline")
  sf_target <- locked %>% filter(.data$cohort == .env$cohort, analysis == "natural_cubic_spline",
    outcome == "recovery_event")
  max_error <- max(abs(coef(sf$fit)[sf_target$term] - sf_target$estimate))
  test <- regTermTest(sf$fit, as.formula(paste("~", paste(nl_names, collapse = "+"))))
  p_target <- nonlinear %>% filter(.data$cohort == .env$cohort, outcome == "recovery")
  stopifnot(max_error < 1e-8, abs(test$p - p_target$p_nonlinear) < 1e-10)

  transform_nl <- function(x) {
    pred_b <- predict(b, newx = x)
    residual_b <- pred_b - cbind(1, x) %*% coef(residual_model)
    residual_b[, q$pivot[seq_len(k)], drop = FALSE] %*%
      solve(qr.R(q)[seq_len(k), seq_len(k), drop = FALSE])
  }
  stopifnot(max(abs(transform_nl(z$delta10) - nl)) < 1e-8)
  grid <- seq(quantile(z$delta10, .01), quantile(z$delta10, .99), length.out = 151)
  beta <- coef(sf$fit)
  variance <- vcov(sf$fit)
  weights <- z$analysis_weight / sum(z$analysis_weight)
  tt <- delete.response(terms(sf$fit))
  rows <- lapply(grid, function(x) {
    nd <- z
    nd$delta10 <- x
    nlx <- drop(transform_nl(x))
    for (j in seq_len(k)) nd[[nl_names[j]]] <- nlx[j]
    mm <- model.matrix(tt, nd)[, names(beta), drop = FALSE]
    p <- plogis(drop(mm %*% beta))
    probability <- sum(weights * p)
    gradient <- colSums(mm * (weights * p * (1 - p)))
    se <- sqrt(drop(t(gradient) %*% variance %*% gradient))
    data.frame(cohort = cohort, signed_difference = x * 10, probability = probability,
      ci_lower = probability - qnorm(.975) * se,
      ci_upper = probability + qnorm(.975) * se)
  })
  curves[[cohort]] <- bind_rows(rows)
  samples[[cohort]] <- data.frame(cohort = cohort, signed_difference = z$delta10 * 10)
  qs <- quantile(z$delta10 * 10, c(.10, .50, .90), names = FALSE)
  quantiles[[cohort]] <- data.frame(cohort = cohort, P10 = qs[1], P50 = qs[2], P90 = qs[3],
    unit = "mL/min/1.73 m2", unit_of_observation = "eligible interval; unweighted")
  flows[[cohort]] <- data.frame(cohort = cohort,
    difficulty_start_intervals = nrow(raw), unknown_next_wave = sum(raw$next_state == "unknown"),
    known_next_wave = sum(raw$next_state != "unknown"),
    excluded_marker_covariate_weight = sum(raw$next_state != "unknown") - nrow(z),
    final_intervals = nrow(z), final_participants = n_distinct(z$pid),
    recovery = sum(z$next_state == "independent"),
    persistent_difficulty = sum(z$next_state == "difficulty"), death = sum(z$next_state == "dead"))
  checks[[cohort]] <- data.frame(cohort = cohort,
    linear_coefficient_absolute_error = abs(coef(f$fit)["delta10"] - target$estimate),
    spline_max_coefficient_absolute_error = max_error,
    nonlinearity_p_absolute_error = abs(test$p - p_target$p_nonlinear),
    spline_basis_reconstruction_max_error = max(abs(transform_nl(z$delta10) - nl)))
}
write.csv(bind_rows(curves), file.path(out, "existing_recovery_spline_curves.csv"), row.names = FALSE)
write.csv(bind_rows(samples), file.path(out, "risk_set_exposure_distribution.csv"), row.names = FALSE)
write.csv(bind_rows(quantiles), file.path(out, "risk_set_exposure_quantiles.csv"), row.names = FALSE)
write.csv(bind_rows(flows), file.path(out, "risk_set_flow.csv"), row.names = FALSE)
write.csv(bind_rows(checks), file.path(out, "existing_model_reproduction_checks.csv"), row.names = FALSE)
write_json(list(
  purpose = "Visualisation of the existing analysis; no new model specification or inferential tests",
  original_analysis = "analysis/code/102_clinical_nonlinear_interaction_competing_risk.R",
  locked_results = "analysis/tables/clinical_enhancement_model_results.csv",
  reproduction = "Recovery linear and spline coefficients, sample sizes, events and nonlinear P values checked against archived outputs",
  curve = "Survey-standardised recovery probabilities at fixed exposure values; existing coefficient covariance and delta method; observed P1-P99 range",
  distribution = "Unweighted eligible intervals; baseline exposure is repeated for participants contributing several intervals"
), file.path(out, "clinical_visual_provenance.json"), pretty = TRUE, auto_unbox = TRUE)
print(bind_rows(flows)); print(bind_rows(quantiles)); print(bind_rows(checks))
