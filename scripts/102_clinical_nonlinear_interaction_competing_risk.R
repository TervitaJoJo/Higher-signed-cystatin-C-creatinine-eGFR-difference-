Sys.setlocale("LC_CTYPE", "English_United States.utf8")
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(survey); library(splines)
  library(metafor); library(jsonlite)
})
options(survey.lonely.psu = "adjust", survey.adjust.domain.lonely = TRUE)
set.seed(20261007)

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) normalizePath(args[1], winslash = "/") else getwd()
dir.create(file.path(root, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(root, "audit"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(root, "reports"), recursive = TRUE, showWarnings = FALSE)

cc <- readRDS(file.path(root, "data_internal", "charls_cohorts.rds"))
hc <- readRDS(file.path(root, "data_internal", "hrs_renal_cohorts.rds"))

canonical_charls <- function(co) {
  p <- co$panel %>% filter(wave <= 4) %>% mutate(
    state = recode(state_difficulty, dependent = "difficulty"),
    psu = factor(community), stratum = factor("CHARLS"),
    draw_lag_months = 0, blood_wt_original = blood_wt,
    date = as.Date(ifelse(is.na(iwy) | is.na(iwm), NA,
      sprintf("%04d-%02d-15", iwy, iwm)))
  )
  list(panel = p, baseline = co$baseline, first_wave = co$first_wave)
}

cohorts <- list(CHARLS2011 = canonical_charls(cc$baseline2011),
                HRS2016VBS = hc$HRS2016VBS)

intervals_for <- function(co) {
  co$panel %>% arrange(pid, wave) %>% group_by(pid) %>% mutate(
    start_state = state, start_severity = difficulty_n,
    next_state = lead(state), next_wave = lead(wave), next_year = lead(year),
    next_date = lead(date), baseline_common_state = first(state),
    prior_recovery = cumsum(as.integer(lag(state, default = "unknown") == "difficulty" &
      state == "independent")) > 0,
    always_independent = cumall(state == "independent")
  ) %>% ungroup() %>% filter(
    baseline_common_state %in% c("independent", "difficulty"),
    !is.na(next_wave), next_wave == wave + 1,
    start_state %in% c("independent", "difficulty")
  ) %>% mutate(
    interval = factor(paste0(year, "-", next_year)),
    age_c = (age - 70) / 10,
    visit_gap_years = as.numeric(next_date - date) / 365.25,
    known = as.integer(next_state != "unknown")
  )
}

make_design <- function(d, is_hrs) {
  if (is_hrs) svydesign(ids = ~psu, strata = ~stratum, weights = ~analysis_weight,
                         data = d, nest = TRUE)
  else svydesign(ids = ~psu, weights = ~analysis_weight, data = d)
}

prepare <- function(co, cohort) {
  is_hrs <- startsWith(cohort, "HRS")
  d <- intervals_for(co) %>% filter(start_state == "difficulty") %>% mutate(
    delta10 = (egfr_cys - egfr_cr) / 10,
    mean_egfr = (egfr_cys + egfr_cr) / 2,
    ckd = as.integer(mean_egfr < 60),
    bmi_c = bmi - 25,
    analysis_weight = blood_wt,
    recovery_event = case_when(next_state == "independent" ~ 1,
      next_state %in% c("difficulty", "dead") ~ 0, TRUE ~ NA_real_),
    persistent_event = case_when(next_state == "difficulty" ~ 1,
      next_state %in% c("independent", "dead") ~ 0, TRUE ~ NA_real_),
    death_event = case_when(next_state == "dead" ~ 1,
      next_state %in% c("independent", "difficulty") ~ 0, TRUE ~ NA_real_)
  ) %>% filter(is.finite(delta10), is.finite(mean_egfr), is.finite(analysis_weight),
               analysis_weight > 0)
  if (!is_hrs) d$interval <- droplevels(d$interval)
  d
}

base_terms <- function(d, cohort) {
  x <- c("delta10", "mean_egfr", "age_c", "female", "education", "interval", "start_severity",
          "bmi", "smoking", "hypertension", "diabetes", "heart", "stroke")
  if (cohort == "HRS2016VBS") x <- c(x, "draw_lag_months")
  x[x %in% names(d)]
}

fit_svy <- function(d, cohort, outcome, terms, label, keep_delta = TRUE) {
  is_hrs <- cohort == "HRS2016VBS"
  needed <- unique(c(outcome, all.vars(reformulate(terms)), "pid", "psu", "analysis_weight", if (is_hrs) "stratum"))
  z <- d[complete.cases(d[, needed, drop = FALSE]), , drop = FALSE]
  z <- z[is.finite(z$analysis_weight) & z$analysis_weight > 0, , drop = FALSE]
  z$analysis_weight <- z$analysis_weight / mean(z$analysis_weight)
  empty <- data.frame(cohort = cohort, analysis = label, outcome = outcome,
    n_people = length(unique(z$pid)), n_intervals = nrow(z), events = sum(z[[outcome]] == 1),
    term = NA_character_, estimate = NA_real_, se = NA_real_, odds_ratio = NA_real_,
    ci_lower = NA_real_, ci_upper = NA_real_, p = NA_real_, note = "insufficient data",
    stringsAsFactors = FALSE)
  if (nrow(z) < 60 || sum(z[[outcome]] == 1, na.rm = TRUE) < 30 ||
      sum(z[[outcome]] == 0, na.rm = TRUE) < 30) return(list(rows = empty, fit = NULL, data = z))
  des <- make_design(z, is_hrs)
  form <- reformulate(terms, response = outcome)
  fit <- tryCatch(svyglm(form, design = des, family = quasibinomial()), error = function(e) e)
  if (inherits(fit, "error")) {
    empty$note <- conditionMessage(fit)
    return(list(rows = empty, fit = NULL, data = z))
  }
  cf <- summary(fit)$coefficients
  use_terms <- if (keep_delta) intersect(rownames(cf), terms) else setdiff(rownames(cf), "(Intercept)")
  use_terms <- use_terms[is.finite(cf[use_terms, 1])]
  rows <- lapply(use_terms, function(term) {
    ci <- suppressMessages(confint(fit, parm = term))
    data.frame(cohort = cohort, analysis = label, outcome = outcome,
      n_people = length(unique(z$pid)), n_intervals = nrow(z),
      events = sum(z[[outcome]] == 1, na.rm = TRUE), term = term,
      estimate = unname(cf[term, 1]), se = unname(cf[term, 2]),
      odds_ratio = exp(unname(cf[term, 1])), ci_lower = exp(ci[1]), ci_upper = exp(ci[2]),
      p = unname(cf[term, 4]), note = "survey-weighted quasibinomial", stringsAsFactors = FALSE)
  })
  list(rows = bind_rows(rows), fit = fit, data = z)
}

model_rows <- list(); model_cache <- list(); spline_rows <- list(); interaction_rows <- list()
outcomes <- c(recovery_event = "recovery", persistent_event = "persistent_difficulty", death_event = "death")
for (cohort in names(cohorts)) {
  d <- prepare(cohorts[[cohort]], cohort)
  bt <- base_terms(d, cohort)
  for (outcome in names(outcomes)) {
    label <- outcomes[[outcome]]
    f <- fit_svy(d, cohort, outcome, bt, "primary_linear")
    model_rows[[length(model_rows) + 1L]] <- f$rows
    model_cache[[paste(cohort, outcome, "linear", sep = "|")]] <- f

    # Natural cubic spline basis with three nonlinear degrees of freedom.
    z <- f$data
    if (!is.null(f$fit) && nrow(z) >= 60) {
      b <- ns(z$delta10, df = 4)
      br <- residuals(lm(b ~ z$delta10))
      q <- qr(br)
      k <- min(3L, q$rank)
      if (k >= 2L) {
        nl <- qr.Q(q)[, seq_len(k), drop = FALSE]
        nl_names <- paste0("rcs_nl", seq_len(k))
        colnames(nl) <- nl_names
        z[, nl_names] <- nl
        st <- c(bt, nl_names)
        sf <- fit_svy(z, cohort, outcome, st, "natural_cubic_spline")
        model_rows[[length(model_rows) + 1L]] <- sf$rows
        nt <- tryCatch(regTermTest(sf$fit, as.formula(paste("~", paste(nl_names, collapse = "+")))),
                       error = function(e) NULL)
        spline_rows[[length(spline_rows) + 1L]] <- data.frame(
          cohort = cohort, outcome = label, n_people = nrow(z),
          nonlinear_df = k, p_nonlinear = if (is.null(nt)) NA_real_ else as.numeric(nt$p),
          p_overall = if (is.null(nt)) NA_real_ else as.numeric(nt$p),
          knots = "natural cubic spline df=4; nonlinear basis orthogonalized to delta10",
          stringsAsFactors = FALSE)
      }
    }
  }

  # Prespecified effect-modification analyses. These are heterogeneity tests, not subgroup screening.
  modifiers <- c(age = "age_c", sex = "female", BMI = "bmi_c", CKD = "ckd",
                 baseline_severity = "start_severity")
  for (mod_name in names(modifiers)) {
    mod <- modifiers[[mod_name]]
    if (!mod %in% names(d)) next
    int_term <- paste0("delta10:", mod)
    terms <- unique(c(bt, mod, int_term))
    for (outcome in names(outcomes)) {
      f <- fit_svy(d, cohort, outcome, terms, paste0("interaction_", mod_name))
      model_rows[[length(model_rows) + 1L]] <- f$rows %>% filter(grepl("delta10", term))
      if (!is.null(f$fit)) {
        cf <- summary(f$fit)$coefficients
        tt <- intersect(c(int_term, paste0(mod, ":delta10")), rownames(cf))
        if (length(tt)) interaction_rows[[length(interaction_rows) + 1L]] <- data.frame(
          cohort = cohort, outcome = outcomes[[outcome]], modifier = mod_name,
          interaction_term = tt[1], estimate = cf[tt[1], 1], se = cf[tt[1], 2],
          odds_ratio = exp(cf[tt[1], 1]), p = cf[tt[1], 4],
          n_people = length(unique(f$data$pid)), n_intervals = nrow(f$data), stringsAsFactors = FALSE)
      }
    }
  }

  # Standardized absolute risk difference for delta10 P90 versus P10.
  for (outcome in names(outcomes)) {
    f <- model_cache[[paste(cohort, outcome, "linear", sep = "|")]]
    if (is.null(f$fit)) next
    z <- f$data; qx <- quantile(z$delta10, c(.10, .90), na.rm = TRUE, names = FALSE)
    tt <- delete.response(terms(f$fit)); x0 <- model.matrix(tt, z)
    zlo <- z; zhi <- z; zlo$delta10 <- qx[1]; zhi$delta10 <- qx[2]
    xlo <- model.matrix(tt, zlo); xhi <- model.matrix(tt, zhi)
    beta <- coef(f$fit); V <- vcov(f$fit)
    keep <- intersect(colnames(xlo), names(beta)); xlo <- xlo[, keep, drop = FALSE]; xhi <- xhi[, keep, drop = FALSE]
    beta <- beta[keep]; V <- V[keep, keep, drop = FALSE]
    pl <- plogis(drop(xlo %*% beta)); ph <- plogis(drop(xhi %*% beta))
    w <- z$analysis_weight / sum(z$analysis_weight)
    risk_low <- sum(w * pl); risk_high <- sum(w * ph); rd <- risk_high - risk_low
    grad <- colSums((ph * (1 - ph) * xhi - pl * (1 - pl) * xlo) * w)
    se <- sqrt(drop(t(grad) %*% V %*% grad)); ci <- rd + c(-1, 1) * qnorm(.975) * se
    beneficial <- if (outcome == "recovery_event") rd > 0 else rd < 0
    measure_label <- ifelse(beneficial, "NNT", "NNH")
    model_rows[[length(model_rows) + 1L]] <- data.frame(
      cohort = cohort, analysis = "standardized_p90_vs_p10", outcome = outcomes[[outcome]],
      n_people = length(unique(z$pid)), n_intervals = nrow(z), events = sum(z[[outcome]] == 1),
      term = "absolute_risk_difference", estimate = rd, se = se, odds_ratio = NA_real_,
      ci_lower = ci[1], ci_upper = ci[2], p = 2 * pnorm(-abs(rd / se)),
      risk_low = risk_low, risk_high = risk_high, nnt_or_nnh = 1 / abs(rd),
      nnt_type = measure_label, exposure_low = qx[1], exposure_high = qx[2],
      note = paste0("delta10 P90=", round(qx[2], 3), " vs P10=", round(qx[1], 3),
                    "; next observed visit risk; ", measure_label, "=", round(1 / abs(rd), 1)),
      stringsAsFactors = FALSE)
  }
}

model_results <- bind_rows(model_rows)
spline_results <- bind_rows(spline_rows)
interaction_results <- bind_rows(interaction_rows)
write.csv(model_results, file.path(root, "tables", "clinical_enhancement_model_results.csv"), row.names = FALSE)
write.csv(spline_results, file.path(root, "tables", "clinical_nonlinear_spline_results.csv"), row.names = FALSE)
write.csv(interaction_results, file.path(root, "tables", "clinical_interaction_results.csv"), row.names = FALSE)

# Random-effects synthesis: two independent cohorts only, hence explicitly exploratory.
meta_rows <- list()
for (outcome_raw in names(outcomes)) {
  outcome <- outcomes[[outcome_raw]]
  z <- model_results %>% filter(analysis == "primary_linear", outcome == !!outcome_raw,
                                term == "delta10", is.finite(estimate), is.finite(se), se > 0)
  if (nrow(z) == 2) {
    m <- tryCatch(rma.uni(yi = z$estimate, sei = z$se, method = "REML"), error = function(e) NULL)
    if (!is.null(m)) {
      pr <- predict(m)
      meta_rows[[length(meta_rows) + 1L]] <- data.frame(
        outcome = outcome, n_cohorts = nrow(z), pooled_or = exp(as.numeric(m$b)),
        ci_lower = exp(m$ci.lb), ci_upper = exp(m$ci.ub), p = m$pval,
        tau2 = m$tau2, I2_percent = m$I2, Q = m$QE, Q_p = m$QEp,
        prediction_lower = exp(pr$pi.lb), prediction_upper = exp(pr$pi.ub),
        scope = "REML random-effects; two independent cohorts; exploratory",
        stringsAsFactors = FALSE)
    }
  }
}
meta_results <- bind_rows(meta_rows)
write.csv(meta_results, file.path(root, "tables", "clinical_random_effects_meta.csv"), row.names = FALSE)

write_json(list(
  status = "completed",
  primary_exposure = "delta10=(eGFRcys-eGFRcr)/10",
  cohorts = names(cohorts),
  nonlinear = "natural cubic spline df=4 with nonlinear basis orthogonalized to linear exposure",
  interactions = names(modifiers),
  competing_risk = "discrete-time cause-specific transition logistic models on difficulty-start intervals",
  absolute_risk = "survey-standardized P90 versus P10 next-visit risk difference; NNT/NNH only descriptive",
  random_effects = "REML with I2 and prediction interval; only two independent cohorts",
  limitation = "visit intervals are discrete and do not provide exact event times; this is not a Fine-Gray event-time analysis"
), file.path(root, "audit", "clinical_enhancement_specification.json"), pretty = TRUE, auto_unbox = TRUE)

report <- c(
  "# Clinical enhancement analysis",
  "",
  "Exposure: delta10=(eGFRcys-eGFRcr)/10. The risk set is an interval beginning with an observed ADL difficulty state; the next scheduled visit is the outcome assessment.",
  "",
  "## Methods",
  "- Natural cubic spline (four degrees of freedom) with a nonlinear basis orthogonalized to delta10; the nonlinear Wald test is reported separately.",
  "- Prespecified interaction tests: age, sex, BMI, CKD defined by mean eGFR <60, and baseline difficulty severity.",
  "- Discrete-time cause-specific transition models for recovery, persistent difficulty, and death. Fine-Gray was not used because exact event dates are unavailable.",
  "- Absolute risks are standardized to the observed covariate distribution for delta10 P10 versus P90 at the next scheduled visit. NNT/NNH is descriptive and should not be interpreted causally.",
  "- Random-effects synthesis is limited to CHARLS2011 and HRS2016VBS; I2 and prediction intervals are exploratory with k=2.",
  "",
  "## Files",
  "- tables/clinical_enhancement_model_results.csv",
  "- tables/clinical_nonlinear_spline_results.csv",
  "- tables/clinical_interaction_results.csv",
  "- tables/clinical_random_effects_meta.csv",
  "- audit/clinical_enhancement_specification.json",
  "",
  "## Interpretation boundary",
  "These additions strengthen dose-response, heterogeneity, competing-outcome and clinical-effect reporting. They do not establish that the muscle oxidative-phosphorylation signal mediates the eGFR discordance association, because the clinical cohorts and muscle transcriptome datasets are not patient-linked."
)
writeLines(report, file.path(root, "reports", "clinical_enhancement_report_zh.md"), useBytes = TRUE)
cat("CLINICAL_ENHANCEMENT_COMPLETE\n")
print(spline_results)
print(interaction_results)
print(meta_results)
