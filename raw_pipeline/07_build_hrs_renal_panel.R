Sys.setlocale("LC_CTYPE", "English_United States.utf8")
suppressPackageStartupMessages({library(haven);library(dplyr);library(tidyr);library(jsonlite)})
args <- commandArgs(trailingOnly=TRUE)
root <- if(length(args)) args[[1]] else getwd()
src <- Sys.getenv("HRS_RAW_DIR")
rand_path <- Sys.getenv("HRS_RAND_FILE")
tracker_path <- Sys.getenv("HRS_TRACKER_FILE")
vbs_path <- Sys.getenv("HRS_VBS_FILE")
dbs_path <- Sys.getenv("HRS_DBS_FILE")
source(file.path(root,"code/renal_equations.R"))
num <- function(x) as.numeric(x)
date_value <- function(x) {
  if(inherits(x,"Date")) return(x)
  fmt <- attr(x,"format.stata")
  if(is.null(fmt)) stop("Interview date lacks a declared Stata date format")
  z <- num(x);ok <- is.finite(z);out <- rep(as.Date(NA),length(z))
  if(grepl("^%tm",fmt)) {
    out[ok] <- as.Date(sprintf("%04d-%02d-15",1960+floor(z[ok]/12),z[ok]%%12+1))
  } else if(grepl("^%td",fmt)) {
    out[ok] <- as.Date(z[ok],origin="1960-01-01")
  } else {
    # RAND codebook PDF p.117: SAS days since 1960-01-01; day is usually imputed as 15.
    out[ok] <- as.Date(z[ok],origin="1960-01-01")
    if(!all(as.integer(format(out[ok],"%Y")) %in% 1990:2025)) stop("RAND documented-date range validation failed")
  }
  out
}
pid9 <- function(hh,pn) {
  a <- as.integer(as.character(hh));b <- as.integer(as.character(pn))
  stopifnot(!anyNA(a),!anyNA(b),all(a>=0 & a<=999999),all(b>=0 & b<=999))
  paste0(sprintf("%06d",a),sprintf("%03d",b))
}
audit <- list()
record <- function(cohort,measure,n) audit[[length(audit)+1L]] <<- data.frame(cohort=cohort,measure=measure,n=n)
items <- c("dress","bath","eat","bed","toilt")
variables <- c("hhidpn","rahhidpn","ragender","raeduc","raestrat","raehsamp","raracem","rahispan","radyear","radmonth",
  unlist(lapply(10:16,function(w) paste0("r",w,c("iwstat","iwmid","agey_m","proxy",paste0(items,"a"),paste0(items,"h"))))),
  unlist(lapply(c(10,13),function(w) paste0("r",w,c("bmi","smoken","hibpe","diabe","hearte","stroke")))))
h <- read_dta(rand_path,col_select=all_of(variables))
h$join_id <- sprintf("%09.0f",num(h$hhidpn))
stopifnot(!anyDuplicated(h$join_id),all(h$join_id==as.character(h$rahhidpn)))
tr <- read_dta(tracker_path,col_select=all_of(c("HHID","PN","PVBSWGTR","MBIOWGTR","STRATUM","SECU"))) %>%
  transmute(join_id=pid9(HHID,PN),vbs_wt_current=num(PVBSWGTR),dbs_wt_current=num(MBIOWGTR),
    tracker_stratum=num(STRATUM),tracker_secu=num(SECU))
stopifnot(!anyDuplicated(tr$join_id))
vbs <- read_sav(vbs_path,col_select=all_of(c("HHID","PN","PCR","PCYSC","PVBSWGTR","PVBS_N_DAYS","PVBSWHY0WGT"))) %>%
  transmute(join_id=pid9(HHID,PN),scr=num(PCR),cys=num(PCYSC),cys_uncalibrated=num(PCYSC),
    blood_wt_original=num(PVBSWGTR),draw_lag_days=num(PVBS_N_DAYS),zero_weight_reason=num(PVBSWHY0WGT))
dbs <- read_dta(dbs_path,col_select=all_of(c("HHID","PN","MCYSC_ADJ","MCYSCUW","MBIOWGTR"))) %>%
  transmute(join_id=pid9(HHID,PN),scr=NA_real_,cys=num(MCYSC_ADJ),cys_uncalibrated=num(MCYSCUW),
    blood_wt_original=num(MBIOWGTR),draw_lag_days=NA_real_,zero_weight_reason=NA_real_)
stopifnot(!anyDuplicated(vbs$join_id),!anyDuplicated(dbs$join_id))
state_from <- function(m) ifelse(rowSums(m==1,na.rm=TRUE)>0,1L,ifelse(rowSums(!is.na(m))==ncol(m),0L,NA_integer_))
id_baselines <- list()
make_cohort <- function(blood,first,last,name,marker,weight_col) {
  record(name,"blood_source_rows",nrow(blood))
  record(name,"blood_matched_rand",sum(blood$join_id %in% h$join_id))
  record(name,"blood_unmatched_rand",sum(!blood$join_id %in% h$join_id))
  linked <- blood %>% inner_join(h,by="join_id",relationship="one-to-one") %>% left_join(tr,by="join_id",relationship="one-to-one")
  prefix <- paste0("r",first)
  record(name,"interview_date_stata_month_format",as.integer(grepl("^%tm",attr(linked[[paste0(prefix,"iwmid")]],"format.stata"))))
  record(name,"interview_date_documented_SAS_days1960",1)
  iw <- date_value(linked[[paste0(prefix,"iwmid")]])
  nextiw <- date_value(linked[[paste0("r",first+1,"iwmid")]])
  age <- num(linked[[paste0(prefix,"agey_m")]])
  female <- case_when(num(linked$ragender)==1 ~ 0, num(linked$ragender)==2 ~ 1,TRUE ~ NA_real_)
  draw_date <- if(first==13) iw+linked$draw_lag_days else iw
  death_end <- rep(as.Date(NA),nrow(linked))
  dy <- num(linked$radyear);dm <- num(linked$radmonth)
  valid_death <- !is.na(dy)&!is.na(dm)&dy>=1990&dy<=2025&dm %in% 1:12
  death_end[valid_death] <- as.Date(sprintf("%04d-%02d-01",dy[valid_death],dm[valid_death]))+31
  # Dates are approximate; only large contradictions warrant exclusion.
  gross_death_conflict <- !is.na(death_end)&!is.na(draw_date)&death_end+31<draw_date
  draw_after_next <- if(first==13) !is.na(nextiw)&!is.na(draw_date)&draw_date>nextiw+31 else rep(FALSE,nrow(linked))
  record(name,"gross_death_before_draw_conflict",sum(gross_death_conflict))
  record(name,"draw_after_first_followup",sum(draw_after_next))
  record(name,"missing_tracker_weight",sum(is.na(linked[[weight_col]])))
  record(name,"positive_current_weight",sum(linked[[weight_col]]>0,na.rm=TRUE))
  record(name,"positive_original_weight",sum(linked$blood_wt_original>0,na.rm=TRUE))
  record(name,"current_original_weight_differ",sum(abs(linked[[weight_col]]-linked$blood_wt_original)>1e-5,na.rm=TRUE))
  record(name,"rand_tracker_stratum_disagreement",sum(num(linked$raestrat)!=linked$tracker_stratum,na.rm=TRUE))
  record(name,"rand_tracker_secu_disagreement",sum(num(linked$raehsamp)!=linked$tracker_secu,na.rm=TRUE))
  base <- data.frame(join_id=linked$join_id,age=age,age_at_blood=age+if(first==13) linked$draw_lag_days/365.25 else 0,
    female=female,education=case_when(num(linked$raeduc) %in% c(1,2) ~ "less_hs",num(linked$raeduc)==3 ~ "hs",
      num(linked$raeduc) %in% c(4,5) ~ "college",TRUE ~ NA_character_),
    bmi_raw=num(linked[[paste0(prefix,"bmi")]]),smoking=num(linked[[paste0(prefix,"smoken")]]),
    hypertension=num(linked[[paste0(prefix,"hibpe")]]),diabetes=num(linked[[paste0(prefix,"diabe")]]),
    heart=num(linked[[paste0(prefix,"hearte")]]),stroke=num(linked[[paste0(prefix,"stroke")]]),
    stratum=num(linked$raestrat),secu=num(linked$raehsamp),
    race=num(linked$raracem),hispanic=num(linked$rahispan),
    scr=linked$scr,cys=linked$cys,cys_uncalibrated=linked$cys_uncalibrated,
    blood_wt=linked[[weight_col]],blood_wt_original=linked$blood_wt_original,
    draw_lag_days=linked$draw_lag_days,draw_lag_months=if(first==13) linked$draw_lag_days/(365.25/12) else 0,
    start_iwstat=num(linked[[paste0(prefix,"iwstat")]]),gross_conflict=gross_death_conflict|draw_after_next,
    draw_date=draw_date)
  for(v in c("smoking","hypertension","diabetes","heart","stroke")) base[[v]][!base[[v]] %in% 0:1] <- NA_real_
  base <- base %>% mutate(education=factor(education,levels=c("less_hs","hs","college")),
    bmi=ifelse(bmi_raw>=12 & bmi_raw<=60,bmi_raw,NA_real_),
    egfr_cr=egfr_cr(scr,age_at_blood,female),egfr_cys=egfr_cys(cys,age_at_blood,female),
    egfr_combined=egfr_combined(scr,cys,age_at_blood,female),
    egfr_cys_uncalibrated=egfr_cys(cys_uncalibrated,age_at_blood,female),
    low_cr=ifelse(is.na(egfr_cr),NA,as.integer(egfr_cr<60)),
    marker_extreme=(!is.na(scr)&(scr<.2|scr>20))|(!is.na(cys)&(cys<.2|cys>10)),
    psu=interaction(stratum,secu,drop=TRUE))
  record(name,"alive_age60_blood_linked",sum(base$age>=60 & base$start_iwstat %in% c(1,2),na.rm=TRUE))
  base <- base %>% filter(age>=60,start_iwstat %in% c(1,2),!gross_conflict,!is.na(.data[[marker]]))
  record(name,"age60_marker_timing_valid",nrow(base))
  p <- bind_rows(lapply(first:last,function(w) {
    mat <- sapply(items,function(j) {
      x <- num(linked[[paste0("r",w,j,"a")]])
      case_when(x==0 ~ 0L,x %in% c(1,2) ~ 1L,TRUE ~ NA_integer_)
    })
    hm <- sapply(items,function(j) {
      x <- num(linked[[paste0("r",w,j,"h")]])
      case_when(x==0 ~ 0L,x %in% c(1,2) ~ 1L,TRUE ~ NA_integer_)
    })
    st <- num(linked[[paste0("r",w,"iwstat")]])
    fun <- state_from(mat)
    data.frame(join_id=linked$join_id,wave=w,year=1990+w*2,iwstat=st,
      date=date_value(linked[[paste0("r",w,"iwmid")]]),proxy=num(linked[[paste0("r",w,"proxy")]]),
      difficulty=fun,difficulty_n=ifelse(rowSums(!is.na(mat))==5,rowSums(mat,na.rm=TRUE),NA_real_),
      help_received=state_from(hm),
      state=case_when(st %in% c(3,5,6) ~ "dead",st %in% c(1,2)&fun==0 ~ "independent",
        st %in% c(1,2)&fun==1 ~ "difficulty",TRUE ~ "unknown"))
  })) %>% arrange(join_id,wave) %>% group_by(join_id) %>%
    mutate(prior_dead=lag(cummax(as.integer(iwstat %in% c(3,5,6))),default=0)) %>% ungroup()
  conflict <- p %>% filter(prior_dead==1,iwstat %in% c(1,2)) %>% distinct(join_id)
  record(name,"alive_after_reported_death",nrow(conflict))
  p$state[p$prior_dead==1] <- "dead"
  first_state <- p %>% filter(wave==first) %>% select(join_id,baseline_state=state,baseline_severity=difficulty_n)
  base <- base %>% left_join(first_state,by="join_id",relationship="one-to-one") %>%
    filter(baseline_state %in% c("independent","difficulty"),!join_id %in% conflict$join_id) %>% mutate(pid=row_number())
  record(name,"included_cohort",nrow(base))
  record(name,"baseline_any_difficulty",sum(base$baseline_state=="difficulty"))
  record(name,"positive_weight_included",sum(base$blood_wt>0,na.rm=TRUE))
  panel <- base %>% inner_join(p,by="join_id",relationship="one-to-many") %>% select(-join_id) %>% arrange(pid,wave)
  id_baselines[[name]] <<- base
  base <- base %>% select(-join_id)
  write.csv(panel %>% count(year,state,name="n"),file.path(root,paste0("tables/",name,"_renal_wave_states.csv")),row.names=FALSE)
  missing <- base %>% summarise(across(c(age,age_at_blood,female,education,bmi,smoking,hypertension,diabetes,heart,stroke,
    stratum,secu,scr,cys,egfr_cr,egfr_cys,egfr_combined,blood_wt,baseline_severity),~sum(is.na(.x)))) %>%
    pivot_longer(everything(),names_to="variable",values_to="n_missing") %>% mutate(n_total=nrow(base),fraction=n_missing/n_total)
  write.csv(missing,file.path(root,paste0("tables/",name,"_renal_missingness.csv")),row.names=FALSE)
  dist <- bind_rows(lapply(c("scr","cys","egfr_cr","egfr_cys","egfr_combined","bmi","blood_wt","draw_lag_days"),function(v) {
    x<-base[[v]][is.finite(base[[v]])]; q<-if(length(x)) unname(quantile(x,c(0,.01,.25,.5,.75,.99,1))) else rep(NA_real_,7)
    data.frame(variable=v,n=length(x),min=q[1],p01=q[2],p25=q[3],median=q[4],p75=q[5],p99=q[6],max=q[7])
  }))
  write.csv(dist,file.path(root,paste0("tables/",name,"_renal_distributions.csv")),row.names=FALSE)
  list(baseline=base,panel=panel,first_wave=first,last_wave=last,primary_marker=marker)
}
co <- list(HRS2016VBS=make_cohort(vbs,13,16,"HRS2016VBS","egfr_cr","vbs_wt_current"),
           HRS2010DBS=make_cohort(dbs,10,14,"HRS2010DBS","egfr_cys","dbs_wt_current"))
common_blood <- intersect(vbs$join_id,dbs$join_id)
record("cross_hrs_samples","blood_source_people_overlap",length(common_blood))
write.csv(bind_rows(audit),file.path(root,"audit/hrs_renal_linkage_and_flow.csv"),row.names=FALSE)
saveRDS(co,file.path(root,"data_internal/hrs_renal_cohorts.rds"))
used <- c(rand_path,tracker_path,vbs_path,dbs_path)
write.csv(data.frame(source=used,md5=unname(tools::md5sum(used)),bytes=file.info(used)$size),file.path(root,"audit/hrs_renal_used_source_hashes.csv"),row.names=FALSE,fileEncoding="UTF-8")
writeLines(capture.output(sessionInfo()),file.path(root,"logs/R_session_hrs_renal_build.txt"))
print(bind_rows(audit))
message("Completed HRS kidney panels. Internal cache has project sequence IDs, not raw person IDs.")

private <- co
for(nm in names(private)) private[[nm]]$baseline <- id_baselines[[nm]]
dir.create(file.path(root,"data_internal/private_linkage"),recursive=TRUE,showWarnings=FALSE)
saveRDS(private,file.path(root,"data_internal/private_linkage/hrs_baseline_keys.rds"))
