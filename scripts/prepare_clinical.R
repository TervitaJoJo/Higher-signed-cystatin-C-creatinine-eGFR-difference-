suppressWarnings(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
suppressPackageStartupMessages({library(dplyr);library(tidyr);library(survey);library(splines);library(jsonlite)})
options(survey.lonely.psu="adjust",survey.adjust.domain.lonely=TRUE)
root <- normalizePath(commandArgs(TRUE)[1],winslash="/")
bundle <- readRDS(file.path(root,"input/clinical_panels.rds"))
saveRDS(list(baseline2011=bundle$CHARLS2011),file.path(root,"data_internal/charls_cohorts.rds"))
saveRDS(list(HRS2016VBS=bundle$HRS2016VBS),file.path(root,"data_internal/hrs_renal_cohorts.rds"))
src <- readLines(file.path(root,"code/102_clinical_nonlinear_interaction_competing_risk.R"))
first <- grep('^cc <- readRDS',src)[1]
last <- grep('^model_rows <- list',src)[1]-1L
eval(parse(text=src[first:last]))
cache <- list(); common <- list()
selection <- readRDS(file.path(root,"input/selection_intervals.rds"))
for (cohort in names(selection)) common[[paste0(cohort,"_egfr_cr")]] <- list(data=selection[[cohort]])
for (cohort in names(cohorts)) {
  d <- prepare(cohorts[[cohort]],cohort)
  f <- fit_svy(d,cohort,"recovery_event",base_terms(d,cohort),"primary_linear")
  z <- f$data
  z$outcome_recovery <- z$recovery_event
  z$outcome_death <- z$death_event
  z$outcome_alive_recovery <- ifelse(z$next_state=="dead",NA,z$recovery_event)
  cache[[paste0(cohort,"_outcome_recovery_continuous")]] <- list(data=z)
  common[[paste0(cohort,"_egfr_cys")]] <- list(data=z)
}
saveRDS(cache,file.path(root,"data_internal/marker_discordance_models.rds"))
saveRDS(common,file.path(root,"data_internal/common_endpoint_models.rds"))
cat("PORTABLE_CLINICAL_INPUTS_PREPARED\n")
