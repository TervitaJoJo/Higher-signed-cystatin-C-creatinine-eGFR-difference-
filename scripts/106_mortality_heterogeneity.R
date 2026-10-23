Sys.setlocale("LC_CTYPE", "English_United States.utf8")
suppressPackageStartupMessages({library(dplyr); library(survey); library(jsonlite)})
options(survey.lonely.psu = "adjust", survey.adjust.domain.lonely = TRUE)
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
root <- normalizePath(args[1], winslash = "/", mustWork = TRUE)
set.seed(20261008)
out <- file.path(root, "mortality_heterogeneity_20261008")
dir.create(file.path(out, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out, "audit"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out, "reports"), recursive = TRUE, showWarnings = FALSE)

cc <- readRDS(file.path(root, "data_internal", "charls_cohorts.rds"))
hc <- readRDS(file.path(root, "data_internal", "hrs_renal_cohorts.rds"))
canonical_charls <- function(co) {
  p <- co$panel %>% filter(wave <= 4) %>% mutate(
    state = recode(state_difficulty, dependent = "difficulty"),
    psu = factor(community), stratum = factor("CHARLS"), draw_lag_months = 0,
    date = as.Date(ifelse(is.na(iwy) | is.na(iwm), NA,
      sprintf("%04d-%02d-15", iwy, iwm)))
  )
  list(panel = p)
}
cohorts <- list(CHARLS2011 = canonical_charls(cc$baseline2011), HRS2016VBS = hc$HRS2016VBS)

intervals_for <- function(co) {
  co$panel %>% arrange(pid, wave) %>% group_by(pid) %>% mutate(
    start_state = state, start_severity = difficulty_n,
    next_state = lead(state), next_wave = lead(wave), next_year = lead(year),
    next_date = lead(date), baseline_common_state = first(state)
  ) %>% ungroup() %>% filter(
    baseline_common_state %in% c("independent", "difficulty"),
    !is.na(next_wave), next_wave == wave + 1,
    start_state == "difficulty"
  ) %>% mutate(
    interval = factor(paste0(year, "-", next_year)),
    age_c = (age - 70) / 10,
    # Use scheduled wave/year spacing. Death intervals often have no observed
    # next-visit date, so date differences would differentially misclassify them.
    visit_gap_years = as.numeric(next_year - year),
    delta10 = (egfr_cys - egfr_cr) / 10,
    mean_egfr = (egfr_cys + egfr_cr) / 2,
    bmi_c = bmi - 25,
    ckd = as.integer(mean_egfr < 60),
    analysis_weight = blood_wt,
    death_event = case_when(next_state == "dead" ~ 1,
      next_state %in% c("independent", "difficulty") ~ 0, TRUE ~ NA_real_)
  ) %>% filter(is.finite(delta10), is.finite(mean_egfr), is.finite(analysis_weight),
               analysis_weight > 0, is.finite(death_event))
}

fit_death <- function(d, cohort, label) {
  is_hrs <- cohort == "HRS2016VBS"
  terms <- c("delta10", "mean_egfr", "age_c", "female", "education", "interval", "start_severity",
             "bmi", "smoking", "hypertension", "diabetes", "heart", "stroke")
  if (is_hrs) terms <- c(terms, "draw_lag_months")
  terms <- terms[terms %in% names(d)]
  required <- unique(c("death_event", "psu", "analysis_weight", terms, if (is_hrs) "stratum"))
  z <- d[complete.cases(d[, required, drop = FALSE]) & is.finite(d$analysis_weight) & d$analysis_weight > 0, , drop = FALSE]
  z$analysis_weight <- z$analysis_weight / mean(z$analysis_weight)
  empty <- data.frame(cohort = cohort, subgroup = label, n_intervals = nrow(z),
    n_people = n_distinct(z$pid), events = sum(z$death_event == 1), odds_ratio = NA_real_,
    ci_low = NA_real_, ci_high = NA_real_, p = NA_real_, status = "insufficient", stringsAsFactors = FALSE)
  if (nrow(z) < 100 || sum(z$death_event == 1) < 20 || sum(z$death_event == 0) < 20) return(empty)
  des <- if (is_hrs) svydesign(ids = ~psu, strata = ~stratum, weights = ~analysis_weight, data = z, nest = TRUE)
    else svydesign(ids = ~psu, weights = ~analysis_weight, data = z)
  fit <- tryCatch(svyglm(reformulate(terms, response = "death_event"), design = des, family = quasibinomial()), error = function(e) e)
  if (inherits(fit, "error")) { empty$status <- conditionMessage(fit); return(empty) }
  cf <- summary(fit)$coefficients["delta10", ]
  ci <- suppressMessages(confint(fit, parm = "delta10"))
  data.frame(cohort = cohort, subgroup = label, n_intervals = nrow(z), n_people = n_distinct(z$pid),
             events = sum(z$death_event == 1), odds_ratio = exp(cf[1]), ci_low = exp(ci[1]),
             ci_high = exp(ci[2]), p = cf[4], status = "ok", stringsAsFactors = FALSE)
}

rows <- list(); interaction_rows <- list(); subgroup_data <- list()
for (cohort in names(cohorts)) {
  d <- intervals_for(cohorts[[cohort]])
  med_gap <- median(d$visit_gap_years[is.finite(d$visit_gap_years)], na.rm = TRUE)
  d$gap_group <- ifelse(is.finite(d$visit_gap_years) & d$visit_gap_years < med_gap, "shorter_than_cohort_median", "at_or_above_cohort_median")
  d$age_group <- ifelse(d$age < 75, "interval_start_age_<75", "interval_start_age_>=75")
  for (g in c("gap_group", "age_group")) {
    for (lv in unique(d[[g]])) {
      z <- d[d[[g]] == lv, , drop = FALSE]
      rows[[length(rows) + 1L]] <- fit_death(z, cohort, paste0(g, "=", lv))
      subgroup_data[[paste(cohort, g, lv, sep = "|")]] <- z
    }
  }
  # Formal effect-modification test on each pre-specified heterogeneity variable.
  is_hrs <- cohort == "HRS2016VBS"
  terms <- c("delta10", "mean_egfr", "age_c", "female", "education", "interval", "start_severity",
             "bmi", "smoking", "hypertension", "diabetes", "heart", "stroke")
  if (is_hrs) terms <- c(terms, "draw_lag_months")
  terms <- terms[terms %in% names(d)]
  for (mod in c("gap_group", "age_group")) {
    z <- d[complete.cases(d[, unique(c("death_event", terms, mod, "psu", "analysis_weight", if (is_hrs) "stratum")), drop = FALSE]), , drop = FALSE]
    z$analysis_weight <- z$analysis_weight / mean(z$analysis_weight)
    if (nrow(z) < 100 || sum(z$death_event == 1) < 20) next
    des <- if (is_hrs) svydesign(ids = ~psu, strata = ~stratum, weights = ~analysis_weight, data = z, nest = TRUE)
      else svydesign(ids = ~psu, weights = ~analysis_weight, data = z)
    int_terms <- c(terms, mod, paste0("delta10:", mod))
    fit <- tryCatch(svyglm(reformulate(int_terms, response = "death_event"), design = des, family = quasibinomial()), error = function(e) NULL)
    if (!is.null(fit)) {
      cf <- summary(fit)$coefficients
      term <- grep("^delta10:", rownames(cf), value = TRUE)
      if (length(term)) interaction_rows[[length(interaction_rows) + 1L]] <- data.frame(
        cohort = cohort, modifier = mod, interaction_term = term[1], estimate = cf[term[1], 1],
        odds_ratio = exp(cf[term[1], 1]), se = cf[term[1], 2], p = cf[term[1], 4],
        n_intervals = nrow(z), events = sum(z$death_event == 1), stringsAsFactors = FALSE)
    }
  }
}
stratified <- bind_rows(rows); interactions <- bind_rows(interaction_rows)
if (nrow(interactions)) interactions$p_bh_four <- p.adjust(interactions$p, "BH")
write.csv(stratified, file.path(out, "tables", "mortality_stratified_results.csv"), row.names = FALSE)
write.csv(interactions, file.path(out, "tables", "mortality_heterogeneity_interactions.csv"), row.names = FALSE)
audit <- list(
  outcome = "next-visit death transition",
  stratifiers = c("visit gap below/at cohort median", "baseline age <75 versus >=75"),
  gap_definition = "cohort-specific median among eligible intervals",
  interaction_adjustment = "BH across four interaction tests",
  interpretation = "Exploratory heterogeneity analysis; death remains a secondary competing transition",
  exact_event_time_limitation = TRUE
)
write_json(audit, file.path(out, "audit", "mortality_heterogeneity_specification.json"), pretty = TRUE, auto_unbox = TRUE)
fmt <- function(x, d = 3) ifelse(is.na(x), "NA", formatC(x, format = "f", digits = d))
report <- c(
  "# Mortality heterogeneity audit",
  "",
  "Death was treated as a secondary competing transition. Stratified models used the same survey-weighted covariate set as the primary transition model and examined visit-gap length and age at the start of the analyzed difficulty interval.",
  "",
  "| Cohort | Subgroup | People | Intervals | Deaths | OR per 10-unit delta10 | 95% CI | P |",
  "|---|---|---:|---:|---:|---:|---:|---:|",
  vapply(seq_len(nrow(stratified)), function(i) paste0("| ", stratified$cohort[i], " | ", stratified$subgroup[i], " | ", stratified$n_people[i], " | ", stratified$n_intervals[i], " | ", stratified$events[i], " | ", fmt(stratified$odds_ratio[i]), " | [", fmt(stratified$ci_low[i]), ", ", fmt(stratified$ci_high[i]), "] | ", fmt(stratified$p[i], 4), " |"), character(1)),
  "",
  "Formal interaction tests:",
  "",
  "| Cohort | Modifier | Interaction OR | P | BH-adjusted P |",
  "|---|---|---:|---:|---:|",
  if (nrow(interactions)) vapply(seq_len(nrow(interactions)), function(i) paste0("| ", interactions$cohort[i], " | ", interactions$modifier[i], " | ", fmt(interactions$odds_ratio[i]), " | ", fmt(interactions$p[i], 4), " | ", fmt(interactions$p_bh_four[i], 4), " |"), character(1)) else "No estimable interaction tests",
  "",
  "These analyses are intended to explain the high mortality heterogeneity, not to promote a mortality subgroup as a new primary finding. Differences in death ascertainment, interval timing, background mortality and cohort composition remain possible explanations and cannot be separated with the available data alone.",
  ""
)
writeLines(report, file.path(out, "reports", "mortality_heterogeneity_zh.md"), useBytes = TRUE)
capture.output(sessionInfo(), file = file.path(out, "audit", "R_sessionInfo.txt"))
cat("MORTALITY_HETEROGENEITY_COMPLETE\n")
print(stratified); print(interactions)
