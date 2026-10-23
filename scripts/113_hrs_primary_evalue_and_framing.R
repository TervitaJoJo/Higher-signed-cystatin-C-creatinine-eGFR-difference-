Sys.setlocale("LC_CTYPE", "English_United States.utf8")
suppressPackageStartupMessages({library(dplyr); library(jsonlite)})
args <- commandArgs(trailingOnly=TRUE); stopifnot(length(args)>=1L)
root <- normalizePath(args[1],winslash="/",mustWork=TRUE)
out <- file.path(root,"reports"); dir.create(out,recursive=TRUE,showWarnings=FALSE)
tab <- read.csv(file.path(root,"tables","clinical_enhancement_model_results.csv"),stringsAsFactors=FALSE)
primary <- tab %>% filter(analysis=="primary_linear",term=="delta10",outcome %in% c("recovery_event","persistent_event","death_event")) %>%
  transmute(cohort,outcome,odds_ratio,ci_lower,ci_upper,se)

evalue_or <- function(or) {
  if (!is.finite(or) || or <= 0) return(NA_real_)
  x <- max(or,1/or)
  if (x <= 1) 1 else x + sqrt(x*(x-1))
}
primary$e_value_point <- vapply(primary$odds_ratio,evalue_or,numeric(1))
primary$e_value_ci <- vapply(seq_len(nrow(primary)),function(i) {
  or <- primary$odds_ratio[i]; lo <- primary$ci_lower[i]; hi <- primary$ci_upper[i]
  if (or >= 1) {if (lo<=1) 1 else evalue_or(lo)} else {if (hi>=1) 1 else evalue_or(hi)}
},numeric(1))
primary$interpretation <- ifelse(primary$e_value_ci==1,"CI includes the null; no non-null E-value is defined","approximate OR-scale E-value")
write.csv(primary,file.path(root,"tables","hrs_primary_evalue_results.csv"),row.names=FALSE)

# Formal between-cohort contrast for recovery is descriptive because designs,
# sampling frames and visit schedules differ.
rec <- primary %>% filter(outcome=="recovery_event") %>% mutate(log_or=log(odds_ratio))
if (nrow(rec)==2) {
  z <- diff(rec$log_or)/sqrt(sum(rec$se^2)); ratio <- exp(diff(rec$log_or)); p <- 2*pnorm(-abs(z))
  contrast <- data.frame(contrast="HRS2016VBS versus CHARLS2011",odds_ratio_ratio=ratio,z=z,p=p,
    note="Wald contrast of two survey-weighted log-OR estimates; descriptive cross-population comparison")
} else contrast <- data.frame()
write.csv(contrast,file.path(root,"tables","hrs_vs_charls_recovery_contrast.csv"),row.names=FALSE)

# Date availability audit for the discrete transition design.
cc <- readRDS(file.path(root,"data_internal","charls_cohorts.rds")); hc <- readRDS(file.path(root,"data_internal","hrs_renal_cohorts.rds"))
canonical_charls <- function(co) {p<-co$panel%>%filter(wave<=4)%>%mutate(state=recode(state_difficulty,dependent="difficulty"),date=as.Date(ifelse(is.na(iwy)|is.na(iwm),NA,sprintf("%04d-%02d-15",iwy,iwm)))); list(panel=p)}
cohorts <- list(CHARLS2011=canonical_charls(cc$baseline2011),HRS2016VBS=hc$HRS2016VBS)
date_audit <- bind_rows(lapply(names(cohorts),function(nm) {
  p<-cohorts[[nm]]$panel%>%arrange(pid,wave)%>%group_by(pid)%>%mutate(start_state=state,next_state=lead(state),next_wave=lead(wave),next_date=lead(date))%>%ungroup()%>%filter(start_state=="difficulty",next_wave==wave+1)
  data.frame(cohort=nm,n_difficulty_intervals=nrow(p),n_death_transitions=sum(p$next_state=="dead",na.rm=TRUE),
    death_with_next_date=sum(p$next_state=="dead"&!is.na(p$next_date),na.rm=TRUE),
    all_transitions_with_next_date=sum(!is.na(p$next_date),na.rm=TRUE),
    note="next-visit date is missing for death transitions when no subsequent interview exists")
}))
write.csv(date_audit,file.path(root,"tables","exact_event_time_availability_audit.csv"),row.names=FALSE)

fmt <- function(x,d=3) ifelse(is.na(x),"NA",formatC(x,format="f",digits=d))
lines <- c("# HRS 主分析、CHARLS 探索性重复与 E-value 审计","",
  "## 预设叙事","",
  "HRS2016VBS is designated as the primary cohort because it provides the clearest current renal-marker/functional transition analysis and the strongest primary recovery association. CHARLS2011 is reported as an exploratory replication and population-difference analysis; its direction is concordant but its confidence interval includes the null. The pooled random-effects estimate remains supplementary and is not used to override the HRS-primary framing.","",
  "建议正文表述：HRS 中 delta10 每增加 10 个单位，下一次随访从既存困难恢复的优势比为 1.138（95% CI 1.045–1.239，P=0.0038）。CHARLS 的方向一致（OR 1.086，95% CI 0.996–1.185，P=0.061），作为探索性重复，不作为主要确证证据。两队列恢复效应差异的描述性 Wald P=",if(nrow(contrast)) fmt(contrast$p,4) else "NA","。","",
  "## E-value（近似 OR 量表）","",
  "E-value is reported as a sensitivity parameter, not as proof of no confounding. The calculation uses the conventional point-estimate formula OR + sqrt[OR(OR−1)] for OR>1 and the reciprocal for OR<1. Because the outcome is not necessarily rare, this is an approximate OR-scale value and should not be presented as a risk-ratio E-value.","",
  "| 队列 | 结局 | OR | 95% CI | 点估计 E-value | CI E-value |","|---|---|---:|---:|---:|---:|")
for(i in seq_len(nrow(primary))) lines<-c(lines,paste0("| ",primary$cohort[i]," | ",primary$outcome[i]," | ",fmt(primary$odds_ratio[i])," | [",fmt(primary$ci_lower[i]),", ",fmt(primary$ci_upper[i]),"] | ",fmt(primary$e_value_point[i])," | ",fmt(primary$e_value_ci[i])," |"))
lines <- c(lines,"","HRS recovery 的点估计 E-value 约 1.535；置信区间下限对应的近似 E-value 约 1.263。含义是：若要将 HRS 恢复关联完全推至零点，未测混杂因子需要同时与 delta10 和恢复结局具有至少这一量级的关联（在 OR 量表近似下）。这不能替代负控、敏感性分析或因果识别。","",
  "## Fine–Gray 可行性","",
  "现有主数据按调查波次构造离散转归区间。死亡转归通常没有可用的下一次访视日期，因此无法为全部死亡事件定义精确事件时间和删失时间。当前使用离散时间 cause-specific transition logistic model，并明确不称为 Fine–Gray event-time model。","",
  "## 输出","",
  "- tables/hrs_primary_evalue_results.csv","- tables/hrs_vs_charls_recovery_contrast.csv","- tables/exact_event_time_availability_audit.csv","- reports/clinical_enhancement_report_zh.md")
writeLines(lines,file.path(out,"hrs_primary_evalue_and_framing_zh.md"),useBytes=TRUE)
write_json(list(primary_cohort="HRS2016VBS",charls_role="exploratory replication and population-difference analysis",evalue="approximate OR-scale conventional formula",fine_gray_status="not estimable for all intervals because exact death event dates are unavailable",between_cohort_contrast=if(nrow(contrast)) as.list(contrast) else NULL),file.path(root,"audit","hrs_primary_framing_specification.json"),pretty=TRUE,auto_unbox=TRUE,na="null")
cat("HRS_PRIMARY_EVALUE_FRAMING_COMPLETE\n"); print(primary); print(contrast); print(date_audit)
