Sys.setlocale("LC_CTYPE", "English_United States.utf8")
suppressPackageStartupMessages({library(haven); library(dplyr); library(survey); library(jsonlite)})
options(survey.lonely.psu = "adjust", survey.adjust.domain.lonely = TRUE)
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
root <- normalizePath(args[1], winslash = "/", mustWork = TRUE)
out <- file.path(root, "renal_extension_20261011")
for (p in c("audit", "tables", "reports", "logs", "data_internal")) {
  dir.create(file.path(out, p), recursive = TRUE, showWarnings = FALSE)
}
source(file.path(root, "code/renal_equations.R"))
cohorts <- c("HRS2016VBS", "CHARLS2011")
outcomes <- c("recovery", "persistent_difficulty", "death")
clinical <- function(cohort) c("age_c", "female", "education", "interval", "start_severity",
  "bmi", "smoking", "hypertension", "diabetes", "heart", "stroke",
  if (cohort == "HRS2016VBS") "draw_lag_months")
design_for <- function(d, cohort) {
  if (cohort == "HRS2016VBS") svydesign(ids = ~psu, strata = ~stratum,
    weights = ~analysis_weight, data = d, nest = TRUE)
  else svydesign(ids = ~psu, weights = ~analysis_weight, data = d)
}
checks <- list(); fits <- list(); coefficient_rows <- list(); diagnostics <- list()
flow <- list(); distributions <- list(); observed <- list(); joint <- list()
risks <- list(); risk_differences <- list(); skipped <- list(); warnings_all <- list()
check <- function(label, value) {
  checks[[length(checks) + 1L]] <<- data.frame(check = label, passed = isTRUE(value))
  if (!isTRUE(value)) stop(label)
}
globalCallingHandlers(warning = function(w) {
  warnings_all[[length(warnings_all) + 1L]] <<- data.frame(message = conditionMessage(w))
  invokeRestart("muffleWarning")
})

write_json(list(date = "2026-10-11", status = "specified_before_new_model_fitting",
  primary_cohort = "HRS2016VBS", exploratory_cohort = "CHARLS2011",
  sample = "Existing complete-case difficulty-start transition risk sets; no new outcome selection",
  filtration_comparison = "delta10 + arithmetic mean eGFR versus delta10 + CKD-EPI 2021 combined eGFR, same sample",
  marker_specs = c("BUN and hemoglobin on their common sample in each cohort",
    "HRS adds serum albumin on the three-marker common sample",
    "BUN uses log2; hemoglobin and albumin use original g/dL units",
    "All positive finite values retained; separate post hoc hemoglobin 5-20 g/dL sensitivity"),
  relative_definition = "ratio=eGFRcys/eGFRcr; low <=0.70, similar >0.70 and <1.30, high >=1.30",
  phenotype = "Combined eGFR >=60, 45-59, <45 crossed with relative discordance; descriptive full grid",
  standardized_recovery = "Low versus similar relative eGFR within combined-eGFR strata; high group excluded; domain survey design; clinical covariates plus combined eGFR; >=100 people and >=20 unique event/nonevent people in each group; no causal interpretation",
  binary_subgroups = "Combined eGFR <60 (combined adjustment), and arithmetic mean <60 (mean adjustment); full-design domain subsetting",
  variance = "HRS supplied PSU/stratum; CHARLS community PSU; original blood weights; repeated intervals nested in PSU",
  exclusions = "No marker-based exclusion from full-sample eGFR comparisons; complete cases only within explicitly labelled marker samples",
  boundaries = c("Single measurement is not confirmed CKD", "No measured-GFR or urine albumin reference",
    "Post hoc exploratory sensitivity analyses", "No mediation or diagnostic accuracy inference",
    "Outcome-specific logistic ORs versus other two outcomes combined; not cause-specific or subdistribution hazards")),
  file.path(out, "audit/analysis_plan.json"), pretty = TRUE, auto_unbox = TRUE)

# Read the exact processed biomarker linkage; raw survey IDs are not required.
cache <- readRDS(file.path(root, "data_internal/marker_discordance_models.rds"))
aug_input <- readRDS(file.path(root, "input/renal_intervals.rds"))
marker_baselines <- lapply(aug_input, function(d) d %>% distinct(pid, .keep_all=TRUE) %>%
  select(pid, bun, log2_bun, hemoglobin, albumin))
source_labels <- data.frame(cohort=c("HRS2016VBS","CHARLS2011"),
  source=c("HRS 2016 VBS processed markers","CHARLS 2011 processed markers"))
fit_one <- function(des, cohort, outcome, predictors, analysis, sample, focal,
                    family_id = NA_character_, group = "all") {
  d <- des$variables
  common <- data.frame(cohort = cohort, outcome = outcome, analysis = analysis, sample = sample, group = group,
    n_people = n_distinct(d$pid), n_intervals = nrow(d), events = sum(d[[outcome]]),
    n_psu = n_distinct(d$psu), n_strata = if (cohort == "HRS2016VBS") n_distinct(d$stratum) else NA_integer_)
  warn <- character()
  f <- tryCatch(withCallingHandlers(svyglm(reformulate(predictors, response = outcome),
    design = des, family = quasibinomial()), warning = function(w) {
      warn <<- c(warn, conditionMessage(w)); invokeRestart("muffleWarning")
    }), error = function(e) e)
  if (inherits(f, "error")) {
    diagnostics[[length(diagnostics) + 1L]] <<- cbind(common, status = "fit_error", note = conditionMessage(f))
    return(NULL)
  }
  mm <- model.matrix(f)
  status <- if (qr(mm)$rank != ncol(mm)) "rank_deficient" else if (!isTRUE(f$converged)) "not_converged" else "ok"
  sd_col <- apply(mm, 2, sd)
  scaled <- scale(mm[, is.finite(sd_col) & sd_col > 1e-10, drop = FALSE])
  diagnostics[[length(diagnostics) + 1L]] <<- cbind(common, status = status, converged = f$converged,
    matrix_rank = qr(mm)$rank, matrix_columns = ncol(mm), df_residual = f$df.residual,
    scaled_condition_number = kappa(scaled, exact = TRUE), max_abs_beta = max(abs(coef(f)), na.rm = TRUE),
    min_fitted = min(fitted(f)), max_fitted = max(fitted(f)), note = paste(unique(warn), collapse = " | "))
  if (status != "ok") return(NULL)
  cf <- summary(f)$coefficients
  for (term in focal) {
    check(paste(cohort, analysis, sample, outcome, group, term, "estimable"), term %in% rownames(cf))
    ci <- suppressMessages(confint(f, parm = term))
    focal_data <- if (startsWith(term, "relative_group")) {
      d[as.character(d$relative_group) == sub("^relative_group", "", term), ]
    } else d
    event_people <- n_distinct(focal_data$pid[focal_data[[outcome]] == 1])
    nonevent_people <- n_distinct(focal_data$pid[focal_data[[outcome]] == 0])
    coefficient_rows[[length(coefficient_rows) + 1L]] <<- cbind(common, term = term,
      n_group_people = n_distinct(focal_data$pid), n_group_intervals = nrow(focal_data),
      group_events = sum(focal_data[[outcome]]), event_people = event_people, nonevent_people = nonevent_people,
      sparse_group = min(event_people, nonevent_people) < 20L,
      log_or = unname(cf[term, 1]), se = unname(cf[term, 2]), odds_ratio = exp(unname(cf[term, 1])),
      ci_lower = exp(ci[1]), ci_upper = exp(ci[2]), p = unname(cf[term, 4]), family = family_id,
      covariates = paste(predictors, collapse = ";"), status = status)
  }
  fits[[paste(cohort, outcome, analysis, sample, group, sep = "__")]] <<- f
  f
}

add_descriptive <- function(d, cohort, scheme) {
  for (unit in c("person", "interval")) {
    z <- if (unit == "person") d %>% arrange(pid, wave) %>% distinct(pid, .keep_all = TRUE) else d
    dz <- design_for(z, cohort)
    groups <- if (scheme == "relative_group") levels(z$relative_group) else
      as.vector(outer(levels(z$egfr_stage), levels(z$relative_group), paste, sep = "__"))
    for (g in groups) {
      ix <- as.character(z[[scheme]]) == g
      a <- z[ix, , drop = FALSE]
      r <- data.frame(cohort = cohort, scheme = scheme, unit = unit, group = g,
        n_people = n_distinct(a$pid), n_observations = nrow(a),
        weighted_share_percent = 100 * sum(a$analysis_weight) / sum(z$analysis_weight),
        n_unique_recovery = n_distinct(a$pid[a$recovery == 1]),
        n_unique_nonrecovery = n_distinct(a$pid[a$recovery == 0]))
      for (y in outcomes) {
        r[[paste0("events_", y)]] <- if (unit == "interval") sum(a[[y]]) else NA_integer_
        r[[paste0("weighted_", y)]] <- NA_real_
        r[[paste0("ci_lower_", y)]] <- NA_real_
        r[[paste0("ci_upper_", y)]] <- NA_real_
        if (unit == "interval" && nrow(a)) {
          sub <- dz[ix, ]
          p <- svymean(reformulate(y), sub)
          r[[paste0("weighted_", y)]] <- unname(coef(p)[1])
          if (length(unique(a[[y]])) == 2L) {
            ci <- suppressMessages(confint(svyciprop(reformulate(y), sub, method = "logit", df = degf(sub))))
            r[[paste0("ci_lower_", y)]] <- ci[1]; r[[paste0("ci_upper_", y)]] <- ci[2]
          }
        }
      }
      observed[[length(observed) + 1L]] <<- r
    }
  }
}

standardize <- function(f, cohort, stage) {
  z <- f$survey.design$variables
  w <- as.numeric(weights(f$survey.design)); w <- w / sum(w)
  beta <- coef(f); V <- vcov(f); crit <- qt(.975, f$df.residual)
  ps <- list(); gs <- list()
  for (g in 0:1) {
    zz <- z; zz$relative_low <- g
    X <- model.matrix(delete.response(terms(f)), zz, contrasts.arg = f$contrasts, xlev = f$xlevels)
    X <- X[, names(beta), drop = FALSE]
    p <- plogis(drop(X %*% beta)); risk <- sum(w * p)
    grad <- colSums(X * (w * p * (1 - p)))
    se <- sqrt(drop(t(grad) %*% V %*% grad))
    ci <- plogis(qlogis(risk) + c(-1, 1) * crit * se / (risk * (1 - risk)))
    risks[[length(risks) + 1L]] <<- data.frame(cohort = cohort, egfr_stage = stage,
      relative_group = ifelse(g == 1, "low", "similar"), target_people = n_distinct(z$pid),
      target_intervals = nrow(z), risk = risk, ci_lower = ci[1], ci_upper = ci[2],
      method = "Survey-coefficient delta method; conditional empirical target distribution")
    ps[[g + 1L]] <- risk; gs[[g + 1L]] <- grad
  }
  diff <- ps[[2]] - ps[[1]]; grad <- gs[[2]] - gs[[1]]
  se <- sqrt(drop(t(grad) %*% V %*% grad))
  risk_differences[[length(risk_differences) + 1L]] <<- data.frame(cohort = cohort, egfr_stage = stage,
    contrast = "low_minus_similar", target_people = n_distinct(z$pid), target_intervals = nrow(z),
    risk_difference = diff, se = se, ci_lower = diff - crit * se, ci_upper = diff + crit * se,
    p = 2 * pt(-abs(diff / se), df = f$df.residual), family = "standardized_recovery_all_cohorts_stages",
    interpretation = "Adjusted next-visit probability contrast, not intervention effect or NNT")
}

expected <- read.csv(file.path(root, "tables/final_competing_transition_results.csv"))
augmented <- list()
for (cohort in cohorts) {
  d <- cache[[paste0(cohort, "_outcome_recovery_continuous")]]$data
  n_initial <- nrow(d)
  d <- d %>% left_join(marker_baselines[[cohort]], by = "pid", relationship = "many-to-one") %>% mutate(
    recovery = as.integer(next_state == "independent"), persistent_difficulty = as.integer(next_state == "difficulty"),
    death = as.integer(next_state == "dead"), delta10 = (egfr_cys - egfr_cr) / 10,
    mean_egfr = (egfr_cys + egfr_cr) / 2, ratio = egfr_cys / egfr_cr,
    relative_group = factor(case_when(ratio <= .70 ~ "low", ratio >= 1.30 ~ "high", TRUE ~ "similar"),
      levels = c("similar", "low", "high")), relative_low = as.integer(relative_group == "low"),
    egfr_stage = factor(case_when(egfr_combined >= 60 ~ "ge60", egfr_combined >= 45 ~ "45to59", TRUE ~ "lt45"),
      levels = c("ge60", "45to59", "lt45")), phenotype = paste(egfr_stage, relative_group, sep = "__"),
    combined_low = as.integer(egfr_combined < 60), mean_low = as.integer(mean_egfr < 60),
    complete_bun_hgb = complete.cases(log2_bun, hemoglobin),
    complete_bun_hgb_alb = complete.cases(log2_bun, hemoglobin, albumin),
    hgb_range_sensitivity = complete_bun_hgb & hemoglobin >= 5 & hemoglobin <= 20,
    analysis_weight = blood_wt / mean(blood_wt))
  check(paste(cohort, "risk set unchanged by raw linkage"), nrow(d) == n_initial)
  check(paste(cohort, "all transitions mutually exclusive and complete"),
    all(d$recovery + d$persistent_difficulty + d$death == 1))
  needed <- c(clinical(cohort), "delta10", "mean_egfr", "egfr_combined", "psu", "stratum", "analysis_weight")
  check(paste(cohort, "primary variables complete"), all(complete.cases(d[, needed])))
  age_measurement <- if (cohort == "HRS2016VBS") d$age_at_blood else d$age
  check(paste(cohort, "combined equation reproduced"),
    max(abs(d$egfr_combined - egfr_combined(d$scr, d$cys, age_measurement, d$female))) < 1e-8)
  check(paste(cohort, "positive finite eGFR ratio"), all(is.finite(d$ratio) & d$ratio > 0))
  check(paste(cohort, "person remains in single PSU"),
    all(d %>% group_by(pid) %>% summarise(n = n_distinct(psu), .groups = "drop") %>% pull(n) == 1L))
  augmented[[cohort]] <- d
  des <- design_for(d, cohort)
  people <- d %>% arrange(pid, wave) %>% distinct(pid, .keep_all = TRUE)
  for (sample in c("full", "bun_hgb", if (cohort == "HRS2016VBS") "bun_hgb_alb", "hgb_5to20")) {
    ix <- switch(sample, full = rep(TRUE, nrow(d)), bun_hgb = d$complete_bun_hgb,
      bun_hgb_alb = d$complete_bun_hgb_alb, hgb_5to20 = d$hgb_range_sensitivity)
    z <- d[ix, ]
    flow[[length(flow) + 1L]] <- data.frame(cohort = cohort, sample = sample,
      n_people = n_distinct(z$pid), n_intervals = nrow(z), excluded_people = n_distinct(d$pid) - n_distinct(z$pid),
      excluded_intervals = nrow(d) - nrow(z), recovery_events = sum(z$recovery),
      persistent_events = sum(z$persistent_difficulty), death_events = sum(z$death))
  }
  for (v in c("bun", "hemoglobin", "albumin", "egfr_combined", "mean_egfr", "ratio")) {
    x <- people[[v]]; valid <- is.finite(x)
    q <- if (any(valid)) quantile(x[valid], c(0, .01, .25, .5, .75, .99, 1), names = FALSE) else rep(NA_real_, 7)
    distributions[[length(distributions) + 1L]] <- data.frame(cohort = cohort, variable = v,
      n_total = nrow(people), n_available = sum(valid), n_missing = sum(!valid),
      min = q[1], p01 = q[2], p25 = q[3], median = q[4], p75 = q[5], p99 = q[6], max = q[7],
      n_hgb_outside_5to20 = if (v == "hemoglobin") sum(valid & (x < 5 | x > 20)) else NA_integer_)
  }
  for (y in outcomes) {
    scope <- ifelse(y == "recovery", "recovery", "secondary")
    f0 <- fit_one(des, cohort, y, c("delta10", "mean_egfr", clinical(cohort)),
      "primary_mean", "full", "delta10")
    check(paste(cohort, y, "primary model available"), !is.null(f0))
    got <- tail(coefficient_rows, 1)[[1]]
    target <- expected[expected$cohort == cohort & expected$outcome ==
      c(recovery = "Recovery", persistent_difficulty = "Persistent difficulty", death = "Death")[[y]], ]
    check(paste(cohort, y, "primary results reproduced"), nrow(target) == 1L &&
      max(abs(c(got$odds_ratio - target$odds_ratio, got$ci_lower - target$ci_lower,
        got$ci_upper - target$ci_upper, got$p - target$p))) < 1e-8)
    fit_one(des, cohort, y, c("delta10", "egfr_combined", clinical(cohort)),
      "primary_combined", "full", "delta10", paste0("combined_adjustment_", scope))
    specifications <- list(
      bun_hgb = list(flag = d$complete_bun_hgb,
        models = list(mean_base = character(), mean_bun = "log2_bun", mean_hgb = "hemoglobin",
          mean_bun_hgb = c("log2_bun", "hemoglobin"), combined_base = character(),
          combined_bun_hgb = c("log2_bun", "hemoglobin"))),
      hgb_5to20 = list(flag = d$hgb_range_sensitivity,
        models = list(mean_base = character(), mean_bun_hgb = c("log2_bun", "hemoglobin"),
          combined_base = character(), combined_bun_hgb = c("log2_bun", "hemoglobin"))))
    if (cohort == "HRS2016VBS") specifications$bun_hgb_alb <- list(flag = d$complete_bun_hgb_alb,
      models = list(mean_base = character(), mean_bun_hgb = c("log2_bun", "hemoglobin"),
        mean_bun_hgb_alb = c("log2_bun", "hemoglobin", "albumin"), combined_base = character(),
        combined_bun_hgb = c("log2_bun", "hemoglobin"), combined_bun_hgb_alb = c("log2_bun", "hemoglobin", "albumin")))
    # Identical two- and three-marker masks share the same baseline and intermediate fits.
    if (cohort == "HRS2016VBS" && identical(d$complete_bun_hgb, d$complete_bun_hgb_alb)) {
      specifications$bun_hgb <- NULL
      specifications$bun_hgb_alb$models$mean_bun <- "log2_bun"
      specifications$bun_hgb_alb$models$mean_hgb <- "hemoglobin"
    }
    for (sample in names(specifications)) {
      s <- specifications[[sample]]; domain <- des[s$flag, ]
      check(paste(cohort, sample, "enough marker-complete people"), n_distinct(domain$variables$pid) > 100)
      for (model in names(s$models)) {
        filtration <- if (startsWith(model, "combined")) "egfr_combined" else "mean_egfr"
        fam <- if (endsWith(model, "base")) NA_character_ else paste0("marker_adjustment_", scope)
        fit_one(domain, cohort, y, c("delta10", filtration, clinical(cohort), s$models[[model]]),
          model, sample, "delta10", fam)
      }
    }
    for (filtration in c("clinical", "mean", "combined")) {
      preds <- c("relative_group", clinical(cohort), switch(filtration,
        clinical = character(), mean = "mean_egfr", combined = "egfr_combined"))
      f <- fit_one(des, cohort, y, preds, paste0("relative3_", filtration), "full",
        c("relative_grouplow", "relative_grouphigh"), paste0("relative_categories_", scope))
      if (!is.null(f)) {
        t <- regTermTest(f, ~relative_group, method = "Wald")
        joint[[length(joint) + 1L]] <- data.frame(cohort = cohort, outcome = y,
          analysis = paste0("relative3_", filtration), test = "relative_group_joint_Wald",
          statistic = as.numeric(t$Ftest), df_num = t$df, df_denom = t$ddf,
          p = as.numeric(t$p), family = paste0("relative_joint_", scope))
      }
    }
    for (definition in c("combined", "mean")) {
      low <- paste0(definition, "_low")
      filtration <- if (definition == "combined") "egfr_combined" else "mean_egfr"
      for (s in 0:1) {
        domain <- des[des$variables[[low]] == s, ]
        fit_one(domain, cohort, y, c("delta10", filtration, clinical(cohort)),
          paste0("binary_", definition, "_domain"), "full", "delta10",
          paste0("binary_domains_", scope), ifelse(s == 0, "ge60", "lt60"))
      }
      f <- fit_one(des, cohort, y, c("delta10", filtration, low, paste0("delta10:", low), clinical(cohort)),
        paste0("binary_", definition, "_interaction"), "full", paste0("delta10:", low),
        "binary_interaction_all")
    }
  }
  add_descriptive(d, cohort, "relative_group")
  add_descriptive(d, cohort, "phenotype")
  for (stage in levels(d$egfr_stage)) {
    ix <- d$egfr_stage == stage & d$relative_group %in% c("similar", "low")
    z <- d[ix, ]
    support <- lapply(c("similar", "low"), function(g) {
      a <- z[z$relative_group == g, ]
      data.frame(group = g, n_people = n_distinct(a$pid),
        event_people = n_distinct(a$pid[a$recovery == 1]), nonevent_people = n_distinct(a$pid[a$recovery == 0]))
    }) %>% bind_rows()
    enough <- n_distinct(z$pid) >= 100 && all(support$event_people >= 20 & support$nonevent_people >= 20)
    if (!enough) {
      skipped[[length(skipped) + 1L]] <- data.frame(cohort = cohort, analysis = "standardized_recovery",
        egfr_stage = stage, n_people = n_distinct(z$pid), n_intervals = nrow(z),
        reason = "Below stated unique-person/event support threshold; descriptive results retained")
      next
    }
    f <- fit_one(des[ix, ], cohort, "recovery", c("relative_low", "egfr_combined", clinical(cohort)),
      "within_stage_relative_recovery", "low_and_similar", "relative_low", "stage_recovery_contrasts", stage)
    if (is.null(f)) {
      skipped[[length(skipped) + 1L]] <- data.frame(cohort = cohort, analysis = "standardized_recovery",
        egfr_stage = stage, n_people = n_distinct(z$pid), n_intervals = nrow(z), reason = "Model not estimable/full rank")
    } else standardize(f, cohort, stage)
  }
}

results <- bind_rows(coefficient_rows)
results$p_bh <- NA_real_
for (fam in na.omit(unique(results$family))) {
  ix <- which(!is.na(results$family) & results$family == fam)
  results$p_bh[ix] <- p.adjust(results$p[ix], "BH")
}
desc <- bind_rows(observed)
joint_tests <- bind_rows(joint) %>% group_by(family) %>% mutate(p_bh = p.adjust(p, "BH")) %>% ungroup()
rd <- bind_rows(risk_differences)
if (nrow(rd)) rd$p_bh <- p.adjust(rd$p, "BH")
changes <- list()
for (cohort in cohorts) for (y in outcomes) {
  for (sample in unique(results$sample[results$cohort == cohort & results$term == "delta10"])) {
    for (filtration in c("mean", "combined")) {
      z <- results %>% filter(.data$cohort == .env$cohort, outcome == y, .data$sample == .env$sample,
        term == "delta10", startsWith(analysis, filtration))
      base <- z[z$analysis == paste0(filtration, "_base"), ]
      if (nrow(base) != 1L) next
      adjusted <- z[!endsWith(z$analysis, "base"), ]
      for (i in seq_len(nrow(adjusted))) {
        a <- adjusted[i, ]
        check(paste(cohort, y, sample, a$analysis, "identical same-sample counts"),
          a$n_people == base$n_people && a$n_intervals == base$n_intervals && a$events == base$events)
        changes[[length(changes) + 1L]] <- data.frame(cohort = cohort, outcome = y, sample = sample,
          filtration = filtration, adjustment = a$analysis, n_people = a$n_people, n_intervals = a$n_intervals,
          base_OR = base$odds_ratio, adjusted_OR = a$odds_ratio, base_logOR = base$log_or,
          adjusted_logOR = a$log_or, logOR_change = a$log_or - base$log_or,
          interpretation = "Descriptive conditional coefficient change; not proportion mediated")
      }
    }
  }
}
diag <- bind_rows(diagnostics)
check("All fitted model diagnostics acceptable", all(diag$status == "ok"))
check("All reported OR confidence intervals valid", all(results$ci_lower <= results$odds_ratio & results$ci_upper >= results$odds_ratio))
check("All reported P values valid", all(is.finite(results$p) & results$p >= 0 & results$p <= 1))
shares <- desc %>% group_by(cohort, scheme, unit) %>% summarise(total = sum(weighted_share_percent), .groups = "drop")
check("Relative and joint phenotype shares sum to 100", all(abs(shares$total - 100) < 1e-8))
check("Both original risk-set sizes preserved", identical(vapply(augmented, nrow, integer(1)), c(HRS2016VBS = 3239L, CHARLS2011 = 2306L)))
if (nrow(rd)) check("Standardized risk difference equals corresponding risks", all(vapply(seq_len(nrow(rd)), function(i) {
  z <- bind_rows(risks) %>% filter(cohort == rd$cohort[i], egfr_stage == rd$egfr_stage[i])
  abs(z$risk[z$relative_group == "low"] - z$risk[z$relative_group == "similar"] - rd$risk_difference[i]) < 1e-10
}, logical(1))))

write.csv(results, file.path(out, "tables/renal_extension_model_results.csv"), row.names = FALSE)
write.csv(bind_rows(flow), file.path(out, "tables/marker_sample_flow.csv"), row.names = FALSE)
write.csv(bind_rows(distributions), file.path(out, "tables/marker_distributions.csv"), row.names = FALSE)
write.csv(bind_rows(changes), file.path(out, "tables/same_sample_coefficient_changes.csv"), row.names = FALSE)
write.csv(desc, file.path(out, "tables/relative_and_joint_phenotype_descriptives.csv"), row.names = FALSE)
write.csv(joint_tests, file.path(out, "tables/relative_discordance_joint_tests.csv"), row.names = FALSE)
write.csv(bind_rows(risks), file.path(out, "tables/standardized_recovery_probabilities.csv"), row.names = FALSE)
write.csv(rd, file.path(out, "tables/standardized_recovery_risk_differences.csv"), row.names = FALSE)
write.csv(bind_rows(skipped), file.path(out, "audit/standardization_support_limits.csv"), row.names = FALSE)
write.csv(diag, file.path(out, "audit/model_diagnostics.csv"), row.names = FALSE)
write.csv(source_labels, file.path(out, "audit/marker_source_labels.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(bind_rows(checks), file.path(out, "audit/validation_checks.csv"), row.names = FALSE)
if (length(warnings_all)) write.csv(bind_rows(warnings_all) %>% count(message, name = "occurrences"),
  file.path(out, "audit/warning_inventory.csv"), row.names = FALSE)
saveRDS(fits, file.path(out, "data_internal/renal_extension_models.rds"))
saveRDS(augmented, file.path(out, "data_internal/renal_extension_intervals.rds"))
writeLines(capture.output(sessionInfo()), file.path(out, "logs/R_session.txt"))
write_json(list(status = "passed", checks = nrow(bind_rows(checks)), model_fits = nrow(diag),
  coefficient_rows = nrow(results), standardized_contrasts = nrow(rd),
  multiplicity = "BH within named coefficient families across cohorts; recovery separated from secondary outcomes; baseline same-sample models descriptive; all sensitivity results retained",
  all_new_analyses_post_hoc = TRUE, raw_identifiers_in_public_outputs = FALSE,
  standardization_uncertainty = "Conditional on empirical target covariates; coefficient survey covariance and residual t degrees of freedom",
  subgroup_variance = "Full original design retained before domain selection",
  manuscript_modified = FALSE), file.path(out, "audit/final_validation.json"), pretty = TRUE, auto_unbox = TRUE)
cat("RENAL_EXTENSION_COMPLETE\n")
print(results %>% filter(outcome == "recovery", analysis %in% c("primary_mean", "primary_combined",
  "mean_bun_hgb", "mean_bun_hgb_alb", "combined_bun_hgb", "combined_bun_hgb_alb", "relative3_combined")) %>%
  select(cohort, analysis, sample, term, n_people, n_intervals, odds_ratio, ci_lower, ci_upper, p, p_bh))
print(rd)
