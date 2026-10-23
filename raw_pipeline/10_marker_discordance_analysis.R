Sys.setlocale("LC_CTYPE", "English_United States.utf8")
suppressPackageStartupMessages({library(dplyr); library(tidyr); library(survey); library(jsonlite)})
options(survey.lonely.psu="adjust", survey.adjust.domain.lonely=TRUE)

args <- commandArgs(trailingOnly=TRUE)
root <- if(length(args)) args[[1]] else getwd()
cache <- readRDS(file.path(root,"data_internal/common_endpoint_models.rds"))
out <- list(); diag <- list(); model_cache <- list()

make_design <- function(d, is_hrs) {
  if(is_hrs) svydesign(ids=~psu, strata=~stratum, weights=~analysis_weight, data=d, nest=TRUE)
  else svydesign(ids=~psu, weights=~analysis_weight, data=d)
}

fit_discordance <- function(d, cohort, exposure, outcome, variant="continuous_delta") {
  is_hrs <- startsWith(cohort,"HRS")
  predictors <- c(exposure,"mean_egfr","age_c","female","education","interval","start_severity",
    "bmi","smoking","hypertension","diabetes","heart","stroke")
  if(cohort=="HRS2016VBS") predictors <- c(predictors,"draw_lag_months")
  predictors <- unique(predictors)
  req <- c(outcome,predictors,"psu","analysis_weight",if(is_hrs) "stratum")
  z <- d[complete.cases(d[,req,drop=FALSE]) & is.finite(d$analysis_weight) & d$analysis_weight>0,,drop=FALSE]
  if(nrow(z)<80 || sum(z[[outcome]]==1)<30 || sum(z[[outcome]]==0)<30) return(NULL)
  z$analysis_weight <- z$analysis_weight/mean(z$analysis_weight)
  des <- make_design(z,is_hrs)
  f <- tryCatch(svyglm(reformulate(predictors,response=outcome), design=des, family=quasibinomial()), error=function(e)e)
  if(inherits(f,"error")) return(NULL)
  cf <- summary(f)$coefficients
  terms <- exposure
  ci <- confint(f,parm=terms)
  b <- unname(cf[terms,"Estimate"]); se <- unname(cf[terms,"Std. Error"])
  exposure_label <- if(exposure=="delta10") "eGFRcys-minus-eGFRcr per 10 units" else "discordance category"
  ans <- data.frame(cohort=cohort, outcome=outcome, variant=variant, exposure=exposure_label,
    n_people=n_distinct(z$pid), n_intervals=nrow(z), events=sum(z[[outcome]]==1),
    mean_egfr_mean=mean(z$mean_egfr), delta_mean=mean(z$delta), delta_sd=sd(z$delta),
    odds_ratio=exp(b), ci_lower=exp(ci[1]), ci_upper=exp(ci[2]), p=cf[terms,"Pr(>|t|)"],
    n_psu=n_distinct(z$psu), n_strata=if(is_hrs) n_distinct(z$stratum) else NA_integer_,
    parameters=length(coef(f)), df_residual=f$df.residual, converged=isTRUE(f$converged),
    matrix_rank=qr(model.matrix(f))$rank, matrix_columns=ncol(model.matrix(f)),
    max_abs_coefficient=max(abs(coef(f))))
  if(exposure=="delta10") ans$interpretation <- "OR per 10-unit higher eGFRcys minus eGFRcr, adjusted for mean eGFR"
  if(exposure=="contrast_cys_only") ans$interpretation <- "eGFRcys <60 and eGFRcr >=60 vs both eGFR >=60; discordant groups only"
  if(exposure=="contrast_cr_only") ans$interpretation <- "eGFRcr <60 and eGFRcys >=60 vs both eGFR >=60; discordant groups only"
  list(result=ans, fit=f, data=z)
}

for(nm in c("CHARLS2011","HRS2016VBS","CHARLS2015")) {
  key <- paste(nm,"egfr_cr",sep="_")
  if(!key %in% names(cache)) next
  d <- cache[[key]]$data
  d <- d %>% mutate(
    delta=egfr_cys-egfr_cr,
    delta10=delta/10,
    mean_egfr=(egfr_cys+egfr_cr)/2,
    age_c=(age-70)/10,
    outcome_recovery=as.numeric(recovery),
    outcome_alive_recovery=ifelse(next_state=="dead",NA_real_,as.numeric(recovery)),
    outcome_death=ifelse(next_state=="unknown",NA_real_,as.numeric(next_state=="dead")),
    discordance_group=case_when(
      egfr_cys>=60 & egfr_cr>=60 ~ "both_high",
      egfr_cys<60 & egfr_cr>=60 ~ "cys_only",
      egfr_cr<60 & egfr_cys>=60 ~ "cr_only",
      egfr_cys<60 & egfr_cr<60 ~ "both_low",
      TRUE ~ NA_character_)) %>%
    filter(is.finite(delta),is.finite(mean_egfr))
  d$discordance_group <- factor(d$discordance_group,levels=c("both_high","cys_only","cr_only","both_low"))
  d$analysis_weight <- d$blood_wt
  # Continuous discordance is the fixed primary analysis; the +/-10 clinically readable cut is secondary.
  for(outcome in c("outcome_recovery","outcome_alive_recovery","outcome_death")) {
    a <- fit_discordance(d,nm,"delta10",outcome,"continuous_delta")
    if(!is.null(a)){out[[length(out)+1L]]<-a$result;model_cache[[paste(nm,outcome,"continuous",sep="_")]]<-a}
    # Threshold sensitivity: each discordant group is compared only with both-high eGFR.
    for(spec in list(list(group="cys_only",name="contrast_cys_only"),list(group="cr_only",name="contrast_cr_only"))) {
      z <- d %>% filter(discordance_group %in% c("both_high",spec$group)) %>%
        mutate(!!spec$name := as.integer(discordance_group==spec$group))
      a <- fit_discordance(z,nm,spec$name,outcome,"threshold_discordance")
      if(!is.null(a)){out[[length(out)+1L]]<-a$result;model_cache[[paste(nm,outcome,spec$name,sep="_")]]<-a}
    }
  }
}

res <- bind_rows(out)
if(nrow(res)) res$p_bh <- p.adjust(res$p,method="BH")
write.csv(res,file.path(root,"tables/marker_discordance_model_results.csv"),row.names=FALSE)
write.csv(res %>% select(cohort, n_people, n_intervals, mean_egfr_mean, delta_mean, delta_sd) %>% distinct(),
  file.path(root,"tables/marker_discordance_sample_summary.csv"),row.names=FALSE)
diag <- res %>% select(cohort,outcome,variant,n_people,n_intervals,events,n_psu,n_strata,parameters,df_residual,converged,matrix_rank,matrix_columns,max_abs_coefficient)
write.csv(diag,file.path(root,"audit/marker_discordance_diagnostics.csv"),row.names=FALSE)
meta <- res %>% filter(variant=="continuous_delta",cohort %in% c("CHARLS2011","HRS2016VBS")) %>%
  mutate(log_or=log(odds_ratio),se=(log(ci_upper)-log(ci_lower))/(2*1.96)) %>%
  group_by(outcome) %>% summarise(n_cohorts=n(),pooled_or=exp(sum(log_or/se^2)/sum(1/se^2)),
    pooled_se=sqrt(1/sum(1/se^2)),ci_lower=exp(sum(log_or/se^2)/sum(1/se^2)-1.96*pooled_se),
    ci_upper=exp(sum(log_or/se^2)/sum(1/se^2)+1.96*pooled_se),
    p=2*pnorm(-abs((sum(log_or/se^2)/sum(1/se^2))/pooled_se)),
    Q=sum((1/se^2)*(log_or-sum(log_or/se^2)/sum(1/se^2))^2),
    I2=ifelse(Q>0,max(0,(Q-(n()-1))/Q)*100,0),.groups="drop")
write.csv(meta,file.path(root,"tables/marker_discordance_meta.csv"),row.names=FALSE)
saveRDS(model_cache,file.path(root,"data_internal/marker_discordance_models.rds"))
write_json(list(
  primary_definition="delta10=(eGFRcys-eGFRcr)/10, continuous, adjusted for mean eGFR and prespecified covariates",
  secondary_definition="eGFRcys <60 and eGFRcr >=60 or the reverse, with concordant >=60 as reference",
  cohorts=unique(res$cohort),
  outcomes=c("recovery including death as no recovery","survivor recovery","death"),
  n_models=nrow(res),
  data_source="Existing complete-case interval cache; no raw source modification"),
  file.path(root,"reports/marker_discordance_analysis_spec.json"),pretty=TRUE,auto_unbox=TRUE)
print(res)
