Sys.setlocale("LC_CTYPE", "English_United States.utf8")
suppressPackageStartupMessages({ library(dplyr); library(survey); library(jsonlite) })
options(survey.lonely.psu = "adjust", survey.adjust.domain.lonely = TRUE)
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 1L)
root <- normalizePath(args[1], winslash = "/", mustWork = TRUE)
out <- file.path(root, "renal_stratification_20261010")
for (part in c("tables", "audit", "reports", "figures", "logs", "data_internal")) {
  dir.create(file.path(out, part), recursive = TRUE, showWarnings = FALSE)
}
input <- file.path(root, "data_internal/marker_discordance_models.rds")
cache <- readRDS(input)
expected <- read.csv(file.path(root, "tables/final_competing_transition_results.csv"))
source(file.path(root, "code/renal_equations.R"))
outcomes <- c("recovery", "persistent_difficulty", "death")
stage_levels <- c("ge60", "45to59", "lt45")
group_levels <- c("both_ge60", "cys_only_lt60", "cr_only_lt60", "both_lt60")
result_rows <- list(); diagnostic_rows <- list(); tests <- list(); model_store <- list()
descriptive_rows <- list(); cross_rows <- list(); flow_rows <- list(); benchmark_rows <- list()
extra_warnings <- list()
globalCallingHandlers(warning = function(w) {
  extra_warnings[[length(extra_warnings) + 1L]] <<- data.frame(
    message = conditionMessage(w), call = paste(deparse(conditionCall(w)), collapse = " "))
  invokeRestart("muffleWarning")
})

clinical_terms <- function(cohort) {
  c("age_c", "female", "education", "interval", "start_severity", "bmi", "smoking",
    "hypertension", "diabetes", "heart", "stroke",
    if (cohort == "HRS2016VBS") "draw_lag_months")
}

design_for <- function(d, cohort) {
  if (cohort == "HRS2016VBS") {
    svydesign(ids = ~psu, strata = ~stratum, weights = ~analysis_weight, data = d, nest = TRUE)
  } else svydesign(ids = ~psu, weights = ~analysis_weight, data = d)
}

describe <- function(d, des, cohort, scheme, unit) {
  levels_all <- if (scheme == "egfr_category") stage_levels else group_levels
  weights <- as.numeric(weights(des))
  rows <- lapply(levels_all, function(g) {
    ix <- as.character(d[[scheme]]) == g
    z <- d[ix, , drop = FALSE]
    w <- weights[ix]
    row <- data.frame(cohort = cohort, scheme = scheme, unit = unit, group = g,
      n_people = n_distinct(z$pid), n_intervals = if (unit == "interval") nrow(z) else NA_integer_,
      n_observations = nrow(z), weighted_share_percent = 100 * sum(w) / sum(weights),
      weighted_mean_egfr_combined = weighted.mean(z$egfr_combined, w),
      weighted_mean_delta = weighted.mean(z$delta, w),
      egfr_combined_min = min(z$egfr_combined), egfr_combined_max = max(z$egfr_combined))
    if (unit == "interval") for (y in outcomes) {
      row[[paste0("events_", y)]] <- sum(z[[y]])
      row[[paste0("people_with_", y)]] <- n_distinct(z$pid[z[[y]] == 1])
      row[[paste0("weighted_", y, "_percent")]] <- 100 * weighted.mean(z[[y]], w)
    }
    row
  })
  bind_rows(rows)
}

fit_model <- function(des, cohort, outcome, terms, analysis, extract_terms, group_var = NULL,
                      fixed_group = "all", unit = "interval") {
  d <- des$variables
  f <- reformulate(terms, response = outcome)
  warnings <- character()
  fit <- tryCatch(withCallingHandlers(
    svyglm(f, design = des, family = quasibinomial()),
    warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") }
  ), error = function(e) e)
  common <- data.frame(cohort = cohort, analysis = analysis, outcome = outcome,
    unit = unit, n_people = n_distinct(d$pid), n_intervals = nrow(d), events = sum(d[[outcome]]),
    n_psu = n_distinct(d$psu), n_strata = if (cohort == "HRS2016VBS") n_distinct(d$stratum) else NA_integer_)
  if (inherits(fit, "error")) {
    diagnostic_rows[[length(diagnostic_rows) + 1L]] <<- cbind(common, group = fixed_group,
      converged = FALSE, status = "fit_error", note = conditionMessage(fit))
    return(NULL)
  }
  mm <- model.matrix(fit)
  col_sd <- apply(mm, 2, sd)
  scaled <- scale(mm[, col_sd > 1e-10, drop = FALSE])
  status <- if (qr(mm)$rank < ncol(mm)) "rank_deficient" else if (!isTRUE(fit$converged)) "not_converged" else "ok"
  diagnostic_rows[[length(diagnostic_rows) + 1L]] <<- cbind(common, group = fixed_group,
    converged = isTRUE(fit$converged), status = status, matrix_rank = qr(mm)$rank,
    matrix_columns = ncol(mm), df_residual = fit$df.residual,
    scaled_condition_number = kappa(scaled, exact = TRUE),
    max_abs_coefficient = max(abs(coef(fit))), min_fitted = min(fitted(fit)),
    max_fitted = max(fitted(fit)), note = paste(unique(warnings), collapse = " | "))
  cf <- summary(fit)$coefficients
  for (term in extract_terms) {
    g <- if (is.null(group_var)) fixed_group else sub(paste0("^", group_var), "", term)
    z <- if (is.null(group_var)) d else d[as.character(d[[group_var]]) == g, , drop = FALSE]
    events_group <- sum(z[[outcome]])
    nonevents_group <- nrow(z) - events_group
    event_people <- n_distinct(z$pid[z[[outcome]] == 1])
    nonevent_people <- n_distinct(z$pid[z[[outcome]] == 0])
    sparse <- min(event_people, nonevent_people) < 20L
    if (!term %in% rownames(cf)) {
      row <- cbind(common, term = term, group = g, n_group_people = n_distinct(z$pid),
        n_group_intervals = nrow(z), group_events = events_group, group_nonevents = nonevents_group,
        people_with_event = event_people, people_without_event = nonevent_people, sparse_group = sparse,
        log_or = NA_real_, se = NA_real_, odds_ratio = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_,
        p = NA_real_, status = "not_estimable")
    } else {
      ci <- suppressMessages(confint(fit, parm = term))
      row <- cbind(common, term = term, group = g, n_group_people = n_distinct(z$pid),
        n_group_intervals = nrow(z), group_events = events_group, group_nonevents = nonevents_group,
        people_with_event = event_people, people_without_event = nonevent_people, sparse_group = sparse,
        log_or = unname(cf[term, 1]), se = unname(cf[term, 2]), odds_ratio = exp(unname(cf[term, 1])),
        ci_lower = exp(ci[1]), ci_upper = exp(ci[2]), p = unname(cf[term, 4]),
        status = if (status != "ok") status else if (sparse) "sparse_exploratory" else "ok")
    }
    result_rows[[length(result_rows) + 1L]] <<- row
  }
  model_store[[paste(cohort, analysis, outcome, fixed_group, sep = "__")]] <<- fit
  fit
}

for (cohort in c("HRS2016VBS", "CHARLS2011")) {
  d <- cache[[paste0(cohort, "_outcome_recovery_continuous")]]$data
  stopifnot(all(d$start_state == "difficulty"), all(d$next_state %in% c("independent", "difficulty", "dead")))
  d <- d %>% mutate(
    recovery = as.integer(next_state == "independent"),
    persistent_difficulty = as.integer(next_state == "difficulty"), death = as.integer(next_state == "dead"),
    delta = egfr_cys - egfr_cr, delta10 = delta / 10, mean_egfr = (egfr_cys + egfr_cr) / 2,
    egfr_category = factor(case_when(egfr_combined >= 60 ~ "ge60", egfr_combined >= 45 ~ "45to59",
      TRUE ~ "lt45"), levels = stage_levels),
    discordance4 = factor(case_when(
      egfr_cr >= 60 & egfr_cys >= 60 ~ "both_ge60", egfr_cr >= 60 & egfr_cys < 60 ~ "cys_only_lt60",
      egfr_cr < 60 & egfr_cys >= 60 ~ "cr_only_lt60", TRUE ~ "both_lt60"), levels = group_levels),
    analysis_weight = blood_wt / mean(blood_wt))
  needed <- c(clinical_terms(cohort), "pid", "psu", "stratum", "analysis_weight", "egfr_combined", "delta10", "mean_egfr")
  stopifnot(all(complete.cases(d[, needed])), all(is.finite(d$egfr_combined)), all(d$egfr_combined > 0),
    all(d$analysis_weight > 0), all(d$recovery + d$persistent_difficulty + d$death == 1))
  age_measurement <- if (cohort == "HRS2016VBS") d$age_at_blood else d$age
  stopifnot(max(abs(d$egfr_combined - egfr_combined(d$scr, d$cys, age_measurement, d$female))) < 1e-8)
  stopifnot(all(d %>% group_by(pid) %>% summarise(n = n_distinct(egfr_combined), .groups = "drop") %>% pull(n) == 1))
  des <- design_for(d, cohort)
  people <- d %>% arrange(pid, wave) %>% distinct(pid, .keep_all = TRUE)
  pdes <- design_for(people, cohort)
  flow_rows[[length(flow_rows) + 1L]] <- data.frame(cohort = cohort, n_people = n_distinct(d$pid),
    n_intervals = nrow(d), recovery_events = sum(d$recovery), persistent_events = sum(d$persistent_difficulty),
    death_events = sum(d$death), n_excluded_from_primary_cache = 0,
    classification_scope = "Fixed blood measurement among people contributing difficulty-start intervals")
  for (scheme in c("egfr_category", "discordance4")) {
    descriptive_rows[[length(descriptive_rows) + 1L]] <- describe(d, des, cohort, scheme, "interval")
    descriptive_rows[[length(descriptive_rows) + 1L]] <- describe(people, pdes, cohort, scheme, "person")
  }
  for (unit in c("person", "interval")) {
    z <- if (unit == "person") people else d
    for (g in group_levels) for (s in stage_levels) {
      ix <- as.character(z$discordance4) == g & as.character(z$egfr_category) == s
      cross_rows[[length(cross_rows) + 1L]] <- data.frame(cohort = cohort, unit = unit,
        discordance_group = g, egfr_category = s, n_people = n_distinct(z$pid[ix]), n_observations = sum(ix),
        weighted_share_percent = 100 * sum(z$analysis_weight[ix]) / sum(z$analysis_weight))
    }
  }
  for (y in outcomes) {
    baseline <- fit_model(des, cohort, y, c("delta10", "mean_egfr", clinical_terms(cohort)),
      "primary_reproduction", "delta10")
    r <- tail(result_rows, 1)[[1]]
    target <- expected[expected$cohort == cohort & expected$outcome ==
      c(recovery = "Recovery", persistent_difficulty = "Persistent difficulty", death = "Death")[[y]], ]
    stopifnot(nrow(target) == 1)
    errors <- abs(c(r$odds_ratio - target$odds_ratio, r$ci_lower - target$ci_lower,
      r$ci_upper - target$ci_upper, r$p - target$p))
    benchmark_rows[[length(benchmark_rows) + 1L]] <- data.frame(cohort = cohort, outcome = y,
      max_absolute_error = max(errors), reproduced = max(errors) < 1e-8)
    stopifnot(max(errors) < 1e-8)

    for (spec in c("discordance4_clinical", "discordance4_mean_adjusted")) {
      terms <- c("discordance4", clinical_terms(cohort), if (spec == "discordance4_mean_adjusted") "mean_egfr")
      f <- fit_model(des, cohort, y, terms, spec, paste0("discordance4", group_levels[-1]), "discordance4")
      tt <- regTermTest(f, ~discordance4, method = "Wald")
      tests[[length(tests) + 1L]] <- data.frame(cohort = cohort, analysis = spec, outcome = y,
        test = "all_four_groups_joint_Wald", df_num = tt$df, df_denom = tt$ddf,
        statistic = as.numeric(tt$Ftest), p = as.numeric(tt$p))
    }
    f <- fit_model(des, cohort, y, c("egfr_category", clinical_terms(cohort)),
      "egfr_category_clinical", paste0("egfr_category", stage_levels[-1]), "egfr_category")
    tt <- regTermTest(f, ~egfr_category, method = "Wald")
    tests[[length(tests) + 1L]] <- data.frame(cohort = cohort, analysis = "egfr_category_clinical", outcome = y,
      test = "all_three_categories_joint_Wald", df_num = tt$df, df_denom = tt$ddf,
      statistic = as.numeric(tt$Ftest), p = as.numeric(tt$p))

    interaction <- fit_model(des, cohort, y,
      c("delta10", "mean_egfr", "egfr_category", "delta10:egfr_category", clinical_terms(cohort)),
      "delta_by_egfr_category_interaction", "delta10")
    tt <- regTermTest(interaction, ~delta10:egfr_category, method = "Wald")
    tests[[length(tests) + 1L]] <- data.frame(cohort = cohort, analysis = "delta_by_egfr_category_interaction",
      outcome = y, test = "two_interaction_terms_joint_Wald", df_num = tt$df, df_denom = tt$ddf,
      statistic = as.numeric(tt$Ftest), p = as.numeric(tt$p))
    for (s in stage_levels) {
      domain <- des[as.character(des$variables$egfr_category) == s, ]
      fit_model(domain, cohort, y, c("delta10", "mean_egfr", clinical_terms(cohort)),
        "delta_within_egfr_category", "delta10", fixed_group = s)
    }
  }
}

results <- bind_rows(result_rows)
results$family <- NA_character_
for (a in unique(results$analysis)) {
  if (a %in% c("primary_reproduction", "delta_by_egfr_category_interaction")) next
  for (scope in c("recovery", "secondary")) {
    ix <- results$analysis == a & if (scope == "recovery") results$outcome == "recovery" else results$outcome != "recovery"
    results$family[ix] <- paste(a, scope, sep = "__")
  }
}
results$p_bh <- NA_real_
for (fam in na.omit(unique(results$family))) {
  ix <- which(!is.na(results$family) & results$family == fam)
  results$p_bh[ix] <- p.adjust(results$p[ix], method = "BH")
}
joint <- bind_rows(tests)
joint$family <- paste0(joint$analysis, "__joint_tests_all_outcomes")
joint <- joint %>% group_by(family) %>% mutate(p_bh = p.adjust(p, method = "BH")) %>% ungroup()
desc <- bind_rows(descriptive_rows); cross <- bind_rows(cross_rows)
diag <- bind_rows(diagnostic_rows); benchmark <- bind_rows(benchmark_rows)
stopifnot(all(benchmark$reproduced), all(diag$converged), all(diag$status == "ok"),
  all(results$p[!is.na(results$p)] >= 0 & results$p[!is.na(results$p)] <= 1),
  all(results$ci_lower <= results$odds_ratio, na.rm = TRUE),
  all(results$odds_ratio <= results$ci_upper, na.rm = TRUE))
share_checks <- desc %>% group_by(cohort, scheme, unit) %>% summarise(total = sum(weighted_share_percent), .groups = "drop")
stopifnot(all(abs(share_checks$total - 100) < 1e-8))
interval_checks <- desc %>% filter(unit == "interval") %>% group_by(cohort, scheme) %>%
  summarise(n = sum(n_intervals), .groups = "drop")
stopifnot(all(interval_checks$n == ifelse(interval_checks$cohort == "HRS2016VBS", 3239, 2306)))

write.csv(results, file.path(out, "tables/renal_model_results.csv"), row.names = FALSE)
write.csv(joint, file.path(out, "tables/renal_joint_tests.csv"), row.names = FALSE)
write.csv(desc, file.path(out, "tables/renal_group_descriptives.csv"), row.names = FALSE)
write.csv(cross, file.path(out, "tables/renal_classification_crosstab.csv"), row.names = FALSE)
write.csv(bind_rows(flow_rows), file.path(out, "tables/renal_sample_flow.csv"), row.names = FALSE)
write.csv(diag, file.path(out, "audit/renal_model_diagnostics.csv"), row.names = FALSE)
write.csv(benchmark, file.path(out, "audit/primary_reproduction.csv"), row.names = FALSE)
write.csv(bind_rows(extra_warnings) %>% count(message, call, name = "occurrences"),
  file.path(out, "audit/variance_warning_inventory.csv"), row.names = FALSE)
saveRDS(model_store, file.path(out, "data_internal/renal_model_cache.rds"))
writeLines(capture.output(sessionInfo()), file.path(out, "logs/R_session.txt"))
write_json(list(
  date = "2026-10-10", status = "completed", primary_cohort = "HRS2016VBS",
  exploratory_cohort = "CHARLS2011", input = input, input_md5 = unname(tools::md5sum(input)),
  code_md5 = unname(tools::md5sum(file.path(root, "code/119_renal_strata_discordance.R"))),
  egfr_equation = "Race-free CKD-EPI 2021 creatinine-cystatin C equation, verified against stored values",
  renal_categories = list(ge60 = "eGFRcr-cys >=60", `45to59` = "45 <= eGFRcr-cys <60", lt45 = "eGFRcr-cys <45"),
  discordance_reference = "both eGFRcr and eGFRcys >=60",
  interpretation_boundary = c("Single measurement; not confirmed CKD or diagnostic reclassification accuracy",
    "No albuminuria or measured-GFR reference in this analysis", "All new analyses are post hoc exploratory",
    "Observed recovery assessed at the next scheduled visit; no exact recovery event time",
    "Separate survey-weighted logistic models compare each transition with the other two combined",
    "These ORs are not Fine-Gray subdistribution or cause-specific hazard ratios",
    "Primary phenotype model adjusts clinical covariates; mean-eGFR adjustment is sensitivity",
    "Do not include delta, mean eGFR, eGFRcr and eGFRcys together: exact linear dependence"),
  survey_variance = "HRS: supplied PSU and stratum; CHARLS: community PSU; repeated intervals nested in PSU",
  domain_variance = "Subset original survey design for eGFR strata, preserving full-design information",
  lonely_psu = "adjust; survey.adjust.domain.lonely=TRUE",
  covariates = clinical_terms("HRS2016VBS"),
  multiple_testing = "BH within each analysis: recovery contrasts across both cohorts; secondary contrasts across both cohorts/outcomes. Joint Wald tests corrected within analysis across all six cohort/outcome tests.",
  sparse_flag = "<20 unique people with events or without events in the focal group; reporting flag, not an exclusion",
  n_model_fits = length(model_store), n_coefficient_rows = nrow(results), n_joint_tests = nrow(joint)
), file.path(out, "audit/analysis_specification.json"), pretty = TRUE, auto_unbox = TRUE)
write_json(list(status = "passed", primary_reproduction = all(benchmark$reproduced),
  group_counts_and_weighted_shares = TRUE, equation_verification = TRUE,
  all_models_full_rank_and_converged = all(diag$status == "ok"), p_values_valid = TRUE,
  confidence_intervals_valid = TRUE, n_fits = length(model_store),
  n_sparse_coefficient_rows = sum(results$sparse_group), patient_rows_in_public_outputs = FALSE),
  file.path(out, "audit/final_validation.json"), pretty = TRUE, auto_unbox = TRUE)
cat("RENAL_STRATIFICATION_COMPLETE\n")
print(desc %>% filter(unit == "person") %>% select(cohort, scheme, group, n_people, weighted_share_percent))
print(results %>% filter(outcome == "recovery", analysis %in% c("discordance4_clinical",
  "discordance4_mean_adjusted", "delta_within_egfr_category")) %>%
  select(cohort, analysis, group, n_group_people, odds_ratio, ci_lower, ci_upper, p, p_bh, status))
print(joint %>% filter(analysis == "delta_by_egfr_category_interaction"))
