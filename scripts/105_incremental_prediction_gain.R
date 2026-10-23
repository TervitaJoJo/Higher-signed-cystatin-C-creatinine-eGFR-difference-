Sys.setlocale("LC_CTYPE", "English_United States.utf8")
suppressPackageStartupMessages({library(dplyr); library(survey); library(jsonlite)})
options(survey.lonely.psu = "adjust", survey.adjust.domain.lonely = TRUE)
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
root <- normalizePath(args[1], winslash = "/", mustWork = TRUE)
set.seed(20261008)
out <- file.path(root, "prediction_gain_20261008")
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
cohorts <- list(CHARLS2011 = canonical_charls(cc$baseline2011),
                HRS2016VBS = hc$HRS2016VBS)

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
    delta10 = (egfr_cys - egfr_cr) / 10,
    mean_egfr = (egfr_cys + egfr_cr) / 2,
    recovery = case_when(next_state == "independent" ~ 1,
      next_state %in% c("difficulty", "dead") ~ 0, TRUE ~ NA_real_),
    analysis_weight = blood_wt
  ) %>% filter(is.finite(delta10), is.finite(mean_egfr),
               is.finite(analysis_weight), analysis_weight > 0)
}

weighted_auc <- function(y, score, w) {
  keep <- is.finite(y) & is.finite(score) & is.finite(w) & w > 0
  y <- y[keep]; score <- score[keep]; w <- w[keep]
  if (sum(y == 1) < 2 || sum(y == 0) < 2) return(NA_real_)
  ord <- order(score, decreasing = FALSE)
  y <- y[ord]; score <- score[ord]; w <- w[ord]
  neg_total <- sum(w[y == 0]); neg_before <- 0; concord <- 0
  i <- 1L
  while (i <= length(y)) {
    j <- i
    while (j < length(y) && score[j + 1L] == score[i]) j <- j + 1L
    wg_pos <- sum(w[i:j][y[i:j] == 1]); wg_neg <- sum(w[i:j][y[i:j] == 0])
    concord <- concord + wg_pos * (neg_before + 0.5 * wg_neg)
    neg_before <- neg_before + wg_neg
    i <- j + 1L
  }
  concord / (sum(w[y == 1]) * neg_total)
}

continuous_nri <- function(y, p0, p1, w) {
  keep <- is.finite(y) & is.finite(p0) & is.finite(p1) & is.finite(w) & w > 0
  y <- y[keep]; p0 <- p0[keep]; p1 <- p1[keep]; w <- w[keep]
  ev <- y == 1; no <- y == 0
  if (!any(ev) || !any(no)) return(NA_real_)
  event_component <- sum(w[ev] * (p1[ev] > p0[ev])) / sum(w[ev]) -
    sum(w[ev] * (p1[ev] < p0[ev])) / sum(w[ev])
  nonevent_component <- sum(w[no] * (p1[no] < p0[no])) / sum(w[no]) -
    sum(w[no] * (p1[no] > p0[no])) / sum(w[no])
  event_component + nonevent_component
}

metric_bootstrap <- function(d, B = 500L) {
  n <- nrow(d); set.seed(20261008 + n)
  vals <- matrix(NA_real_, nrow = B, ncol = 3L)
  for (b in seq_len(B)) {
    ii <- sample.int(n, n, replace = TRUE)
    z <- d[ii, , drop = FALSE]
    vals[b, 1] <- weighted_auc(z$recovery, z$pred_base, z$analysis_weight)
    vals[b, 2] <- weighted_auc(z$recovery, z$pred_plus, z$analysis_weight) - vals[b, 1]
    vals[b, 3] <- continuous_nri(z$recovery, z$pred_base, z$pred_plus, z$analysis_weight)
  }
  data.frame(
    base_auc_low = quantile(vals[, 1], .025, na.rm = TRUE),
    base_auc_high = quantile(vals[, 1], .975, na.rm = TRUE),
    delta_auc_low = quantile(vals[, 2], .025, na.rm = TRUE),
    delta_auc_high = quantile(vals[, 2], .975, na.rm = TRUE),
    nri_low = quantile(vals[, 3], .025, na.rm = TRUE),
    nri_high = quantile(vals[, 3], .975, na.rm = TRUE),
    bootstrap_replicates = B
  )
}

rows <- list(); prediction_cache <- list()
for (cohort_name in names(cohorts)) {
  d <- intervals_for(cohorts[[cohort_name]]) %>% group_by(pid) %>% arrange(wave, .by_group = TRUE) %>% slice(1L) %>% ungroup()
  is_hrs <- cohort_name == "HRS2016VBS"
  base_terms <- c("mean_egfr", "age_c", "female", "education", "interval", "start_severity",
                  "bmi", "smoking", "hypertension", "diabetes", "heart", "stroke")
  if (is_hrs) base_terms <- c(base_terms, "draw_lag_months")
  base_terms <- base_terms[base_terms %in% names(d)]
  required <- c("recovery", "pid", "psu", "analysis_weight", base_terms, if (is_hrs) "stratum")
  d <- d[complete.cases(d[, required, drop = FALSE]) & is.finite(d$analysis_weight) & d$analysis_weight > 0, , drop = FALSE]
  d$analysis_weight <- d$analysis_weight / mean(d$analysis_weight)
  ids <- unique(d$pid)
  set.seed(20261008 + nchar(cohort_name))
  fold_map <- data.frame(pid = ids, fold = sample(rep(1:5, length.out = length(ids))))
  d <- left_join(d, fold_map, by = "pid")
  d$pred_base <- d$pred_plus <- NA_real_
  for (fold in 1:5) {
    tr <- d$fold != fold; te <- d$fold == fold
    if (sum(te) < 20 || sum(d$recovery[tr] == 1) < 20 || sum(d$recovery[tr] == 0) < 20) next
    train <- d[tr, , drop = FALSE]; test <- d[te, , drop = FALSE]
    des <- if (is_hrs) svydesign(ids = ~psu, strata = ~stratum, weights = ~analysis_weight, data = train, nest = TRUE)
      else svydesign(ids = ~psu, weights = ~analysis_weight, data = train)
    f0 <- svyglm(reformulate(base_terms, response = "recovery"), design = des, family = quasibinomial())
    f1 <- svyglm(reformulate(c(base_terms, "delta10"), response = "recovery"), design = des, family = quasibinomial())
    d$pred_base[te] <- as.numeric(predict(f0, newdata = test, type = "response"))
    d$pred_plus[te] <- as.numeric(predict(f1, newdata = test, type = "response"))
  }
  d <- d[is.finite(d$pred_base) & is.finite(d$pred_plus), , drop = FALSE]
  auc0 <- weighted_auc(d$recovery, d$pred_base, d$analysis_weight)
  auc1 <- weighted_auc(d$recovery, d$pred_plus, d$analysis_weight)
  nri <- continuous_nri(d$recovery, d$pred_base, d$pred_plus, d$analysis_weight)
  boot <- metric_bootstrap(d)
  rows[[length(rows) + 1L]] <- data.frame(
    cohort = cohort_name, n_people = nrow(d), events = sum(d$recovery == 1),
    base_auc = auc0, incremental_auc = auc1, delta_auc = auc1 - auc0,
    continuous_nri = nri, base_auc_low = boot$base_auc_low, base_auc_high = boot$base_auc_high,
    delta_auc_low = boot$delta_auc_low, delta_auc_high = boot$delta_auc_high,
    nri_low = boot$nri_low, nri_high = boot$nri_high,
    bootstrap_replicates = boot$bootstrap_replicates,
    base_model = paste(base_terms, collapse = "+"),
    incremental_model = paste(c(base_terms, "delta10"), collapse = "+"),
    note = "5-fold person-level cross-validated predictions; weighted AUC and continuous NRI are exploratory",
    stringsAsFactors = FALSE)
  prediction_cache[[cohort_name]] <- d
}
results <- bind_rows(rows)
write.csv(results, file.path(out, "tables", "incremental_prediction_gain.csv"), row.names = FALSE)
saveRDS(prediction_cache, file.path(out, "data_internal_prediction_cache.rds"))
spec <- list(
  exposure = "delta10=(eGFRcys-eGFRcr)/10",
  base_model = "mean eGFR plus primary clinical covariates",
  incremental_model = "base model plus delta10",
  collinearity_warning = "Do not include eGFRcr, eGFRcys and mean eGFR together; mean eGFR is an exact linear combination of eGFRcr and eGFRcys.",
  outcome = "first observed difficulty interval next-visit recovery; death coded as non-recovery",
  validation = "5-fold person-level cross-validation, survey-weighted training models",
  metrics = c("weighted AUC", "AUC difference", "continuous NRI"),
  limitation = "Exploratory discrimination analysis; not a clinically validated prediction model"
)
write_json(spec, file.path(out, "audit", "incremental_prediction_gain_specification.json"), pretty = TRUE, auto_unbox = TRUE)
report <- c(
  "# Incremental prediction gain audit",
  "",
  "The prediction question was defined as whether adding delta10 improves next-visit recovery discrimination beyond mean eGFR and the prespecified clinical covariates.",
  "",
  "The base model uses mean eGFR. A model containing eGFRcr, eGFRcys and mean eGFR simultaneously is not interpretable because mean eGFR is exactly (eGFRcr + eGFRcys)/2.",
  "",
  "Predictions were generated by five-fold cross-validation at the person level, using each participant's first eligible difficulty interval. Death was coded as non-recovery for this binary discrimination analysis.",
  "",
  "| Cohort | People | Events | Base AUC | Incremental AUC | Delta AUC | Delta AUC bootstrap 95% CI | Continuous NRI (95% CI) |",
  "|---|---:|---:|---:|---:|---:|---:|---:|",
  vapply(seq_len(nrow(results)), function(i) paste0("| ", results$cohort[i], " | ", results$n_people[i], " | ", results$events[i], " | ", round(results$base_auc[i], 3), " | ", round(results$incremental_auc[i], 3), " | ", round(results$delta_auc[i], 3), " | [", round(results$delta_auc_low[i], 3), ", ", round(results$delta_auc_high[i], 3), "] | ", round(results$continuous_nri[i], 3), " [", round(results$nri_low[i], 3), ", ", round(results$nri_high[i], 3), "] |"), character(1)),
  "",
  "These metrics are exploratory and should be described as incremental discrimination, not as proof that delta10 is ready for clinical deployment. The comparison is meaningful relative to mean eGFR; it is not meaningful if the base model already contains both eGFRcr and eGFRcys.",
  ""
)
writeLines(report, file.path(out, "reports", "incremental_prediction_gain_zh.md"), useBytes = TRUE)
capture.output(sessionInfo(), file = file.path(out, "audit", "R_sessionInfo.txt"))
cat("INCREMENTAL_PREDICTION_GAIN_COMPLETE\n")
print(results)
