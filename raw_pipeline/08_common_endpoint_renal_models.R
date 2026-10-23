Sys.setlocale("LC_CTYPE", "English_United States.utf8")
suppressPackageStartupMessages({library(dplyr);library(tidyr);library(survey);library(ggplot2);library(jsonlite)})
options(survey.lonely.psu="adjust",survey.adjust.domain.lonely=TRUE)
set.seed(20261006)
args <- commandArgs(trailingOnly=TRUE)
root <- if(length(args)) args[[1]] else getwd()
hc <- readRDS(file.path(root,"data_internal/hrs_renal_cohorts.rds"))
cc <- readRDS(file.path(root,"data_internal/charls_cohorts.rds"))
canonical_charls <- function(co) {
  p <- co$panel %>% filter(wave<=4) %>% mutate(state=recode(state_difficulty,dependent="difficulty"),
    psu=factor(community),stratum=1,draw_lag_months=0,blood_wt_original=blood_wt,
    date=as.Date(ifelse(is.na(iwy)|is.na(iwm),NA,sprintf("%04d-%02d-15",iwy,iwm))))
  list(panel=p,baseline=co$baseline,first_wave=co$first_wave,last_wave=4,primary_marker="egfr_cr")
}
cohorts <- c(list(CHARLS2011=canonical_charls(cc$baseline2011),CHARLS2015=canonical_charls(cc$landmark2015)),hc)
intervals_for <- function(co) {
  p <- co$panel %>% arrange(pid,wave) %>% group_by(pid) %>% mutate(
    start_state=state,start_severity=difficulty_n,next_state=lead(state),next_wave=lead(wave),next_year=lead(year),
    next_date=lead(date),baseline_common_state=first(state),
    prior_recovery=cumsum(as.integer(lag(state,default="unknown")=="difficulty" & state=="independent"))>0,
    always_independent=cumall(state=="independent")) %>% ungroup() %>%
    filter(baseline_common_state %in% c("independent","difficulty"),!is.na(next_wave),next_wave==wave+1,
      start_state %in% c("independent","difficulty")) %>%
    mutate(interval=factor(paste0(year,"-",next_year)),age_c=(age-70)/10,
      visit_gap_years=as.numeric(next_date-date)/365.25,known=as.integer(next_state!="unknown"),
      recovery=case_when(start_state=="difficulty" & next_state=="independent" ~ 1,
        start_state=="difficulty" & next_state %in% c("difficulty","dead") ~ 0,TRUE ~ NA_real_),
      alive_recovery=ifelse(next_state=="dead",NA,recovery),
      death=ifelse(next_state=="unknown",NA,as.integer(next_state=="dead")),
      onset=case_when(start_state=="independent" & next_state=="difficulty" ~ 1,
        start_state=="independent" & next_state %in% c("independent","dead") ~ 0,TRUE ~ NA_real_))
  p
}
results <- list();diagnostics <- list();counts <- list();flows <- list();ipw_diagnostics <- list();fit_cache <- list()
make_design <- function(d,is_hrs) {
  if(is_hrs) svydesign(ids=~psu,strata=~stratum,weights=~analysis_weight,data=d,nest=TRUE)
  else svydesign(ids=~psu,weights=~analysis_weight,data=d)
}
fit_model <- function(d,cohort,outcome="recovery",marker="egfr_cr",adjustment="full",weighted=TRUE,variant="main",extra_weight=NULL) {
  is_hrs <- startsWith(cohort,"HRS")
  d <- d %>% mutate(exposure=if(marker=="low_cr") low_cr else (90-.data[[marker]])/10)
  predictors <- c("exposure","age_c","female","education")
  if(n_distinct(d$interval)>1) predictors <- c(predictors,"interval")
  if(outcome %in% c("recovery","alive_recovery","death")) predictors <- c(predictors,"start_severity")
  if(adjustment=="full") predictors <- c(predictors,"bmi","smoking","hypertension","diabetes","heart","stroke")
  if(cohort=="HRS2016VBS") predictors <- c(predictors,"draw_lag_months")
  if(variant=="race_hispanic_adjusted") {
    d$race <- factor(d$race)
    predictors <- c(predictors,"race","hispanic")
  }
  if(variant=="bmi_missing_indicator") {
    d$bmi_missing <- as.integer(is.na(d$bmi))
    d$bmi[is.na(d$bmi)] <- median(d$bmi,na.rm=TRUE)
    if(n_distinct(d$bmi_missing)>1) predictors <- c(predictors,"bmi_missing")
  }
  d$analysis_weight <- if(weighted) {
    if(variant=="original_blood_weight") d$blood_wt_original else d$blood_wt
  } else 1
  if(!is.null(extra_weight)) d$analysis_weight <- d$analysis_weight*extra_weight
  n_risk <- nrow(d);unknown <- sum(is.na(d[[outcome]]));n_known <- sum(!is.na(d[[outcome]]))
  required <- c(outcome,predictors,"psu","analysis_weight",if(is_hrs) "stratum")
  d <- d[complete.cases(d[,required]) & is.finite(d$analysis_weight)&d$analysis_weight>0,,drop=FALSE]
  d$interval <- droplevels(d$interval);d$education <- droplevels(d$education)
  events <- sum(d[[outcome]]==1)
  common <- data.frame(cohort=cohort,endpoint="five_ADL_any_difficulty",outcome=outcome,marker=marker,
    adjustment=adjustment,weighted=weighted,variant=variant,n_risk_intervals=n_risk,n_outcome_unknown=unknown,
    n_intervals=nrow(d),n_people=n_distinct(d$pid),events=events,n_known_excluded=n_known-nrow(d),
    n_psu=n_distinct(d$psu),n_strata=if(is_hrs) n_distinct(d$stratum) else NA_integer_,
    n_egfr_below60=if(marker=="low_cr") sum(d$exposure==1) else sum(d[[marker]]<60),
    events_egfr_below60=if(marker=="low_cr") sum(d$exposure==1&d[[outcome]]==1) else sum(d[[marker]]<60&d[[outcome]]==1))
  if(nrow(d)<60 || min(events,nrow(d)-events)<30 || sd(d$exposure)==0) {
    results[[length(results)+1L]] <<- cbind(common,log_or=NA,se=NA,odds_ratio=NA,ci_lower=NA,ci_upper=NA,p=NA,
      note="Insufficient events/variation; descriptive analysis only")
    return(invisible(NULL))
  }
  d$analysis_weight <- d$analysis_weight/mean(d$analysis_weight)
  des <- make_design(d,is_hrs)
  captured <- character()
  fit <- tryCatch(withCallingHandlers(svyglm(reformulate(predictors,response=outcome),design=des,family=quasibinomial()),
    warning=function(w){captured <<- c(captured,conditionMessage(w));invokeRestart("muffleWarning")}),error=function(e)e)
  if(inherits(fit,"error")) {
    results[[length(results)+1L]] <<- cbind(common,log_or=NA,se=NA,odds_ratio=NA,ci_lower=NA,ci_upper=NA,p=NA,note=conditionMessage(fit))
    return(invisible(NULL))
  }
  cf <- summary(fit)$coefficients["exposure",];ci <- confint(fit,parm="exposure")
  results[[length(results)+1L]] <<- cbind(common,log_or=unname(cf[1]),se=unname(cf[2]),odds_ratio=exp(unname(cf[1])),
    ci_lower=exp(ci[1]),ci_upper=exp(ci[2]),p=unname(cf[4]),note=paste(unique(captured),collapse=" | "))
  mm <- model.matrix(fit)
  diagnostics[[length(diagnostics)+1L]] <<- cbind(common,converged=fit$converged,df_residual=fit$df.residual,
    parameters=length(coef(fit)),minority_outcomes_per_parameter=min(events,nrow(d)-events)/length(coef(fit)),
    matrix_rank=qr(mm)$rank,matrix_columns=ncol(mm),max_abs_coefficient=max(abs(coef(fit))),
    min_fitted=min(fitted(fit)),max_fitted=max(fitted(fit)),weight_max=max(d$analysis_weight),
    effective_n_kish=sum(d$analysis_weight)^2/sum(d$analysis_weight^2))
  if(outcome=="recovery" && adjustment=="full" && weighted && variant=="main") {
    fit_cache[[paste(cohort,marker,sep="_")]] <<- list(fit=fit,data=d,cohort=cohort,marker=marker)
  }
  invisible(fit)
}

for(nm in names(cohorts)) {
  co <- cohorts[[nm]];d <- intervals_for(co)
  counts[[length(counts)+1L]] <- d %>% count(interval,start_state,next_state,name="n") %>% mutate(cohort=nm)
  rec <- d %>% filter(start_state=="difficulty")
  primary_marker <- co$primary_marker
  markers <- if(nm=="HRS2010DBS") "egfr_cys" else c("egfr_cr","egfr_cys","egfr_combined")
  for(ex in markers) for(wt in c(FALSE,TRUE)) fit_model(rec,nm,marker=ex,weighted=wt)
  for(wt in c(FALSE,TRUE)) fit_model(rec,nm,marker=primary_marker,adjustment="core",weighted=wt)
  if(nm!="HRS2010DBS") for(wt in c(FALSE,TRUE)) fit_model(rec,nm,marker="low_cr",weighted=wt)
  for(out in c("alive_recovery","death")) fit_model(rec,nm,outcome=out,marker=primary_marker)
  # Added after the first estimates: characterize death and survivor recovery for all marker definitions.
  for(ex in setdiff(markers,primary_marker)) for(out in c("alive_recovery","death")) {
    fit_model(rec,nm,outcome=out,marker=ex,variant="component_decomposition_exploratory")
  }
  onset <- d %>% filter(start_state=="independent",always_independent)
  recurrence <- d %>% filter(start_state=="independent",prior_recovery)
  fit_model(onset,nm,outcome="onset",marker=primary_marker,variant="first_observed_difficulty")
  fit_model(recurrence,nm,outcome="onset",marker=primary_marker,variant="recurrence_after_recovery")
  if(n_distinct(rec$interval)>1) fit_model(rec %>% filter(wave==co$first_wave),nm,marker=primary_marker,variant="first_interval_only")
  fit_model(rec,nm,marker=primary_marker,variant="bmi_missing_indicator")
  if(startsWith(nm,"HRS")) fit_model(rec,nm,marker=primary_marker,variant="original_blood_weight")
  if(nm=="HRS2016VBS") {
    for(ex in markers) {
      fit_model(rec %>% filter(draw_lag_days<=90),nm,marker=ex,variant="draw_within90_days")
      fit_model(rec,nm,marker=ex,variant="race_hispanic_adjusted")
    }
  }
  if(nm=="HRS2010DBS") fit_model(rec,nm,marker="egfr_cys_uncalibrated",variant="DBS_uncalibrated_measurement")
  gap <- d %>% filter(next_state %in% c("independent","difficulty")) %>% group_by(interval) %>%
    summarise(n=n(),n_valid_gap=sum(!is.na(visit_gap_years)),median_years=median(visit_gap_years,na.rm=TRUE),
      min_years=min(visit_gap_years,na.rm=TRUE),max_years=max(visit_gap_years,na.rm=TRUE),.groups="drop")
  write.csv(gap,file.path(root,paste0("audit/",nm,"_common_visit_gaps.csv")),row.names=FALSE)
  raw <- rec %>% count(next_state,name="n") %>% mutate(cohort=nm)
  flows[[length(flows)+1L]] <- raw
  # Sensitivity to observed next-visit availability; these weights do not impute covariates.
  if(startsWith(nm,"HRS")) {
    ap <- c("age_c","female","education","start_severity","bmi","smoking","hypertension","diabetes","heart","stroke")
    if(n_distinct(rec$interval)>1) ap <- c(ap,"interval")
    if(nm=="HRS2016VBS") ap <- c(ap,"draw_lag_months")
    rec$exposure <- (90-rec[[primary_marker]])/10
    required <- c(ap,"exposure","psu","stratum","blood_wt")
    ra <- rec[complete.cases(rec[,required])&rec$blood_wt>0,,drop=FALSE]
    if(sum(ra$known==0)>20) {
      ra$analysis_weight <- ra$blood_wt/mean(ra$blood_wt)
      des <- make_design(ra,TRUE)
      af <- svyglm(reformulate(c("exposure",ap),response="known"),design=des,family=quasibinomial())
      prob <- as.numeric(predict(af,newdata=ra,type="response"))
      iw <- 1/pmax(prob,.05);cap <- unname(quantile(iw,.99));iw <- pmin(iw,cap)
      ipw_diagnostics[[length(ipw_diagnostics)+1L]] <- data.frame(cohort=nm,n_intervals=nrow(ra),n_unknown=sum(ra$known==0),
        min_availability=min(prob),median_availability=median(prob),ipw_cap=cap,n_probability_floor=sum(prob<.05))
      fit_model(ra,nm,marker=primary_marker,variant="availability_IPW_p99",extra_weight=iw)
    }
  }
}
res <- bind_rows(results)
primary_flag <- res$cohort=="HRS2016VBS" & res$outcome=="recovery" & res$marker=="egfr_cr" &
  res$adjustment=="full" & res$weighted & res$variant=="main"
res$p_bh_secondary <- NA_real_
res$p_bh_secondary[!primary_flag] <- p.adjust(res$p[!primary_flag],method="BH")
write.csv(res,file.path(root,"tables/common_endpoint_renal_model_results.csv"),row.names=FALSE)
write.csv(bind_rows(diagnostics),file.path(root,"audit/common_endpoint_model_diagnostics.csv"),row.names=FALSE)
write.csv(bind_rows(ipw_diagnostics),file.path(root,"audit/hrs_renal_availability_weight_diagnostics.csv"),row.names=FALSE)
write.csv(bind_rows(counts),file.path(root,"tables/common_endpoint_transition_counts.csv"),row.names=FALSE)
write.csv(bind_rows(flows),file.path(root,"tables/common_endpoint_recovery_outcome_totals.csv"),row.names=FALSE)
key <- res %>% filter(outcome=="recovery",adjustment=="full",weighted,variant=="main",marker!="low_cr")
write.csv(key,file.path(root,"tables/common_endpoint_key_results.csv"),row.names=FALSE)
saveRDS(fit_cache,file.path(root,"data_internal/common_endpoint_models.rds"))

# Two independent cohorts only. CHARLS2015 and HRS2010 are additional overlapping analyses.
metas <- bind_rows(lapply(c("egfr_cr","egfr_cys","egfr_combined"),function(ex) {
  z <- key %>% filter(cohort %in% c("CHARLS2011","HRS2016VBS"),marker==ex)
  stopifnot(nrow(z)==2,all(is.finite(z$log_or)),all(z$se>0))
  w <- 1/z$se^2;m <- sum(w*z$log_or)/sum(w);se <- sqrt(1/sum(w))
  Q <- sum(w*(z$log_or-m)^2);I2 <- if(Q>0) max(0,(Q-1)/Q)*100 else 0
  data.frame(marker=ex,n_cohorts=2,pooled_or=exp(m),ci_lower=exp(m-1.96*se),ci_upper=exp(m+1.96*se),
    p=2*pnorm(-abs(m/se)),Q=Q,Q_p=pchisq(Q,1,lower.tail=FALSE),I2_percent=I2,
    scope="Exploratory fixed-effect summary; two cohorts with different visit schedules")
}))
metas$p_bh <- p.adjust(metas$p,method="BH")
write.csv(metas,file.path(root,"tables/common_endpoint_exploratory_meta.csv"),row.names=FALSE)
g <- ggplot(key,aes(y=cohort,x=odds_ratio,xmin=ci_lower,xmax=ci_upper,color=cohort))+
  geom_vline(xintercept=1,linetype=2,color="#858A90")+geom_errorbar(orientation="y",width=.18)+geom_point(size=2.8)+
  facet_wrap(~marker,nrow=1)+scale_x_log10()+
  scale_color_manual(values=c(CHARLS2011="#326E8A",CHARLS2015="#279D85",HRS2016VBS="#CE8834",HRS2010DBS="#646D79"))+
  labs(x="Recovery OR per 10 lower eGFR (95% CI)",y=NULL,title="Common endpoint: five ADL difficulties; weighted complete-case models")+
  theme_minimal(base_size=10)+theme(legend.position="none")
ggsave(file.path(root,"figures/common_endpoint_recovery_estimates.png"),g,width=12,height=4.8,dpi=180)
tr <- bind_rows(counts) %>% filter(start_state=="difficulty") %>% group_by(cohort,interval) %>%
  mutate(total=sum(n),percent=100*n/total) %>% ungroup()
tr$interval <- factor(as.character(tr$interval),levels=sort(unique(as.character(tr$interval))))
g2 <- ggplot(tr,aes(x=interval,y=percent,fill=factor(next_state,levels=c("independent","difficulty","dead","unknown"))))+
  geom_col(width=.65)+facet_wrap(~cohort,scales="free_x",ncol=2)+
  scale_fill_manual(values=c(independent="#279D85",difficulty="#DDA638",dead="#535A64",unknown="#BFC6CE"),name="Next visit")+
  labs(x=NULL,y="Percentage of difficulty-start intervals",title="Unweighted observed recovery, persistent difficulty, death and unknown")+
  theme_minimal(base_size=10)+theme(legend.position="bottom",panel.grid.major.x=element_blank())
ggsave(file.path(root,"figures/common_endpoint_observed_transitions.png"),g2,width=11,height=7,dpi=180)
writeLines(capture.output(sessionInfo()),file.path(root,"logs/R_session_common_endpoint.txt"))
write_json(list(hrs_primary=as.list(res[primary_flag,]),n_models=nrow(res),
  endpoint="Five common ADL any-difficulty items; death as no recovery; unknown separately reported",
  independent_cohorts_for_meta=c("CHARLS2011","HRS2016VBS")),file.path(root,"reports/hrs_extension_summary.json"),pretty=TRUE,auto_unbox=TRUE,na="null")
print(key %>% select(cohort,marker,n_people,n_intervals,events,odds_ratio,ci_lower,ci_upper,p))
print(metas)
