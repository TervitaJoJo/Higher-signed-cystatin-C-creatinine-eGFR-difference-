Sys.setlocale("LC_CTYPE", "English_United States.utf8")
suppressPackageStartupMessages({
  library(haven); library(dplyr); library(tidyr); library(jsonlite)
})
args <- commandArgs(trailingOnly=TRUE)
root <- if (length(args)) args[[1]] else getwd()
src <- Sys.getenv("CHARLS_RAW_DIR")
dir.create(file.path(root, "data_internal"), showWarnings=FALSE)
num <- function(x) as.numeric(x)
audit <- list()
record <- function(stage, measure, n) {
  audit[[length(audit)+1L]] <<- data.frame(stage=stage, measure=measure, n=n)
}
used <- character()
read_selected <- function(path, required, optional=character()) {
  full <- file.path(src, path)
  available <- names(read_dta(full, n_max=0))
  stopifnot(all(required %in% available))
  used <<- unique(c(used, full))
  read_dta(full, col_select=all_of(c(required, intersect(optional, available))))
}
unique_ids <- function(d, col="ID", stage) {
  stopifnot(!anyNA(d[[col]]), !any(d[[col]]==""), !anyDuplicated(d[[col]]))
  record(stage, "source_rows", nrow(d))
}

items <- c("dress", "bath", "eat", "bed", "toilt")
covars <- c("agey", "mbmi", "smoken", "hibpe", "diabe", "hearte", "stroke")
h <- read_selected("Harmonized_CHARLS/H_CHARLS_D_Data.dta",
  c("ID", "ID_w1", "communityID", "ragender", "raeduc_c", "radyear", "radmonth",
    unlist(lapply(1:4, function(w) paste0("r", w, c("iwstat", "iwy", "iwm", "agey", paste0(items, "a"))))),
    unlist(lapply(c(1,3), function(w) paste0("r", w, covars[-1])))))
unique_ids(h, stage="harmonized")
idmap <- h %>% select(ID, ID_w1)
map1 <- idmap %>% filter(!is.na(ID_w1), ID_w1!="")
unique_ids(map1, col="ID_w1", stage="wave1_id_map")
link_ids <- function(d, stage) {
  raw <- as.character(d$ID)
  by_current <- match(raw, h$ID)
  by_old <- match(raw, map1$ID_w1)
  linked <- h$ID[by_current]
  fill <- is.na(linked) & !is.na(by_old)
  linked[fill] <- map1$ID[by_old[fill]]
  record(stage, "matched_current_id", sum(!is.na(by_current)))
  record(stage, "matched_wave1_id_only", sum(fill))
  record(stage, "unmatched_id", sum(is.na(linked)))
  d$ID <- linked
  d <- d %>% filter(!is.na(ID))
  unique_ids(d, stage=paste0(stage, "_linked"))
  d
}

source(file.path(root,"code/renal_equations.R"))

b1 <- read_selected("2011/Blood_20140429.dta", c("ID", "newcrea", "cystatinc", "bloodweight"))
unique_ids(b1, stage="blood2011")
b1 <- link_ids(b1, "blood2011") %>% transmute(ID, scr=num(newcrea), cys=num(cystatinc), blood_wt=num(bloodweight))
b3 <- read_selected("2015/Blood.dta", c("ID", "bl_crea", "bl_cysc", "Blood_weight"))
unique_ids(b3, stage="blood2015")
b3 <- link_ids(b3, "blood2015") %>% transmute(ID, scr=num(bl_crea), cys=num(bl_cysc), blood_wt=num(Blood_weight))

base_for <- function(w, blood) {
  d <- h %>% transmute(ID, community=as.character(communityID),
    age=num(.data[[paste0("r",w,"agey")]]),
    female=case_when(num(ragender)==1 ~ 0, num(ragender)==2 ~ 1, TRUE ~ NA_real_),
    # Three prespecified categories avoid an almost empty tertiary-education cell.
    education=case_when(num(raeduc_c)==1 ~ "none", num(raeduc_c) %in% 2:4 ~ "primary_or_less",
                        num(raeduc_c) %in% 5:10 ~ "middle_or_more", TRUE ~ NA_character_),
    bmi=num(.data[[paste0("r",w,"mbmi")]]), smoking=num(.data[[paste0("r",w,"smoken")]]),
    hypertension=num(.data[[paste0("r",w,"hibpe")]]), diabetes=num(.data[[paste0("r",w,"diabe")]]),
    heart=num(.data[[paste0("r",w,"hearte")]]), stroke=num(.data[[paste0("r",w,"stroke")]]),
    alive_start=num(.data[[paste0("r",w,"iwstat")]])==1) %>% left_join(blood, by="ID", relationship="one-to-one")
  for (v in c("smoking", "hypertension", "diabetes", "heart", "stroke")) d[[v]][!d[[v]] %in% 0:1] <- NA_real_
  record(paste0("baseline",w),"nonpositive_bmi_invalid",sum(!is.na(d$bmi) & d$bmi<=0))
  d$bmi[!is.na(d$bmi) & d$bmi<=0] <- NA_real_
  d <- d %>% mutate(egfr_cr=egfr_cr(scr,age,female), egfr_cys=egfr_cys(cys,age,female),
                    egfr_combined=egfr_combined(scr,cys,age,female),
                    low_cr=ifelse(is.na(egfr_cr),NA,as.integer(egfr_cr<60)),
                    marker_extreme=(!is.na(scr) & (scr<0.2 | scr>20)) | (!is.na(cys) & (cys<0.2 | cys>10)),
                    bmi_extreme=!is.na(bmi) & (bmi<12 | bmi>60),bmi_raw=bmi,
                    education=factor(education, levels=c("none","primary_or_less","middle_or_more")))
  record(paste0("baseline",w),"alive_age60",sum(d$alive_start & d$age>=60,na.rm=TRUE))
  d <- d %>% filter(alive_start, age>=60)
  record(paste0("baseline",w),"age60_bmi_outside12_60",sum(d$bmi_extreme))
  d$bmi[d$bmi_extreme] <- NA_real_
  record(paste0("baseline",w),"age60_egfr_cr",sum(!is.na(d$egfr_cr)))
  record(paste0("baseline",w),"age60_egfr_cys",sum(!is.na(d$egfr_cys)))
  record(paste0("baseline",w),"age60_egfr_combined",sum(!is.na(d$egfr_combined)))
  record(paste0("baseline",w),"marker_extreme_flag",sum(d$marker_extreme))
  d
}
base1 <- base_for(1,b1)
base3 <- base_for(3,b3)

state_from <- function(mat) {
  positive <- rowSums(mat==1, na.rm=TRUE)>0
  complete <- rowSums(!is.na(mat))==ncol(mat)
  ifelse(positive, 1L, ifelse(complete, 0L, NA_integer_))
}
count_from <- function(mat) ifelse(rowSums(!is.na(mat))==ncol(mat), rowSums(mat,na.rm=TRUE), NA_real_)
stopifnot(identical(state_from(rbind(c(0,0,0),c(0,NA,0),c(1,NA,0))),c(0L,NA_integer_,1L)))
panel <- list()
raw_audit <- list()
years <- c(2011,2013,2015,2018,2020)
for (w in 1:4) {
  yr <- years[w]
  path <- paste0(yr, if(w==1) "/health_status_and_functioning.dta" else "/Health_Status_and_Functioning.dta")
  raw <- read_selected(path, c("ID", paste0("db0",10:14)), c("proxy","db034"))
  unique_ids(raw, stage=paste0("health",yr))
  raw <- link_ids(raw,paste0("health",yr))
  d <- h %>% select(ID, all_of(paste0("r",w,c("iwstat","iwy","iwm","agey",paste0(items,"a"))))) %>%
    left_join(raw, by="ID", relationship="one-to-one")
  dep <- anydiff <- matrix(NA_real_, nrow(d),5, dimnames=list(NULL,items))
  for (j in 1:5) {
    rv <- num(d[[paste0("db0",9+j)]])
    ha <- num(d[[paste0("r",w,items[j],"a")]])
    valid <- rv %in% 1:4
    # Only a documented harmonized negative can resolve a raw structural skip.
    structural_negative <- !valid & !is.na(ha) & ha==0
    dep[valid,j] <- as.integer(rv[valid]>=3)
    dep[structural_negative,j] <- 0
    anydiff[ha %in% 0:1,j] <- ha[ha %in% 0:1]
    raw_audit[[length(raw_audit)+1L]] <- data.frame(wave=w,year=yr,item=items[j],
      raw_valid=sum(valid),structural_negative_fill=sum(structural_negative),
      raw_harmonized_conflict=sum(valid & !is.na(ha) & ha!=as.integer(rv>=2)),
      dependency_missing=sum(is.na(dep[,j])))
  }
  proxy <- rep(NA_integer_, nrow(d))
  if (w==1) proxy <- ifelse(num(d$proxy) %in% 0:1,num(d$proxy),NA)
  if (w==3) proxy <- case_when(num(d$proxy)==1 ~ 1L, num(d$proxy)==2 ~ 0L, TRUE ~ NA_integer_)
  if ("db034" %in% names(d)) proxy[num(d$db034) %in% 1:4] <- 1L
  panel[[w]] <- data.frame(ID=d$ID,wave=w,year=yr,iwstat=num(d[[paste0("r",w,"iwstat")]]),
    iwy=num(d[[paste0("r",w,"iwy")]]),iwm=num(d[[paste0("r",w,"iwm")]]),
    age_visit=num(d[[paste0("r",w,"agey")]]),proxy=proxy,
    dependency=state_from(dep),severity=count_from(dep),
    difficulty=state_from(anydiff),difficulty_n=count_from(anydiff),
    dependency4=state_from(dep[,-4,drop=FALSE]),severity4=count_from(dep[,-4,drop=FALSE]))
}

s5 <- read_selected("2020/Sample_Infor.dta", c("ID","died","iyear","imonth"))
unique_ids(s5, stage="sample2020")
s5 <- link_ids(s5,"sample2020")
x5 <- read_selected("2020/Exit_Module.dta",c("ID","exb001_1","exb001_2"))
unique_ids(x5,stage="exit2020")
x5 <- link_ids(x5,"exit2020")
d5 <- read_selected("2020/Health_Status_and_Functioning.dta",c("ID",paste0("db00",c(1,3,5,7,9))),c("proxy_5"))
unique_ids(d5,stage="health2020")
d5 <- link_ids(d5,"health2020")
p5 <- h %>% select(ID) %>% left_join(s5,by="ID",relationship="one-to-one") %>%
  left_join(d5,by="ID",relationship="one-to-one") %>% left_join(x5,by="ID",relationship="one-to-one")
dep5 <- diff5 <- matrix(NA_real_,nrow(p5),5)
for (j in 1:5) {
  rv <- num(p5[[paste0("db00",2*j-1)]])
  dep5[rv %in% 1:4,j] <- as.integer(rv[rv %in% 1:4]>=3)
  diff5[rv %in% 1:4,j] <- as.integer(rv[rv %in% 1:4]>=2)
  raw_audit[[length(raw_audit)+1L]] <- data.frame(wave=5,year=2020,item=items[j],raw_valid=sum(rv %in% 1:4),
    structural_negative_fill=0, raw_harmonized_conflict=NA_integer_, dependency_missing=sum(is.na(dep5[,j])))
}
# The 2020 extension uses observed raw responses only; no unverified negative fills.
panel[[5]] <- data.frame(ID=p5$ID,wave=5,year=2020,
  iwstat=case_when(num(p5$died)==1 ~ 5, rowSums(!is.na(dep5))>0 ~ 1, num(p5$died)==0 ~ 4, TRUE ~ 9),
  iwy=num(p5$iyear),iwm=num(p5$imonth),age_visit=NA_real_,
  proxy=case_when(num(p5$proxy_5)==1 ~ 1L, num(p5$proxy_5)==2 ~ 0L, TRUE ~ NA_integer_),
  dependency=state_from(dep5),severity=count_from(dep5),difficulty=state_from(diff5),difficulty_n=count_from(diff5),
  dependency4=state_from(dep5[,-4,drop=FALSE]),severity4=count_from(dep5[,-4,drop=FALSE]))
record("health2020","reported_dead_with_raw_function",sum(num(p5$died)==1 & rowSums(!is.na(dep5))>0,na.rm=TRUE))

long <- bind_rows(panel) %>% left_join(h %>% transmute(ID,death_year=num(radyear),death_month=num(radmonth)),by="ID")
for (endpoint in c("dependency","difficulty","dependency4")) {
  long[[paste0("state_",endpoint)]] <- case_when(long$iwstat %in% c(5,6) ~ "dead",
    long$iwstat==1 & long[[endpoint]]==0 ~ "independent", long$iwstat==1 & long[[endpoint]]==1 ~ "dependent",
    TRUE ~ "unknown")
}
# Death is absorbing. Contradictions are counted and stop the build if present.
absorbing_conflicts <- long %>% arrange(ID,wave) %>% group_by(ID) %>%
  mutate(prior_death=lag(cummax(as.integer(iwstat %in% c(5,6))),default=0)) %>%
  summarise(conflict=any(prior_death==1 & iwstat==1),.groups="drop")
record("panel","alive_after_reported_death",sum(absorbing_conflicts$conflict))
if(any(absorbing_conflicts$conflict)) {
  conflict_ids <- absorbing_conflicts$ID[absorbing_conflicts$conflict]
  # Quarantine conflicting histories instead of silently changing the death date.
  long <- long %>% filter(!ID %in% conflict_ids)
  base1 <- base1 %>% filter(!ID %in% conflict_ids)
  base3 <- base3 %>% filter(!ID %in% conflict_ids)
}
long <- long %>% arrange(ID,wave) %>% group_by(ID) %>%
  mutate(ever_dead=cummax(as.integer(iwstat %in% c(5,6)))) %>% ungroup()
record("panel","prior_death_carried_to_unlisted_wave",sum(long$ever_dead==1 & !long$iwstat %in% c(5,6)))
for (endpoint in c("dependency","difficulty","dependency4")) {
  long[[paste0("state_",endpoint)]][long$ever_dead==1] <- "dead"
}

id_baselines <- list()
make_cohort <- function(base, first_wave) {
  base <- base %>% filter(!is.na(egfr_cr))
  start <- long %>% filter(wave==first_wave) %>% select(ID, starts_with("state_"),severity)
  record(paste0("cohort",first_wave),"before_function_filter",nrow(base))
  base <- base %>% inner_join(start,by="ID",relationship="one-to-one") %>% filter(state_dependency!="unknown",state_dependency!="dead")
  record(paste0("cohort",first_wave),"included",nrow(base))
  record(paste0("cohort",first_wave),"baseline_dependent",sum(base$state_dependency=="dependent"))
  record(paste0("cohort",first_wave),"low_egfr_cr",sum(base$low_cr==1))
  base <- base %>% mutate(pid=row_number())
  allrows <- base %>% select(-starts_with("state_"),-severity) %>% inner_join(long %>% filter(wave>=first_wave),by="ID",relationship="one-to-many")
  allrows <- allrows %>% select(-ID) %>% arrange(pid,wave)
  id_baselines[[as.character(first_wave)]] <<- base
  base <- base %>% select(-ID)
  list(baseline=base,panel=allrows,first_wave=first_wave)
}
cohorts <- list(baseline2011=make_cohort(base1,1),landmark2015=make_cohort(base3,3))
for (nm in names(cohorts)) {
  co <- cohorts[[nm]]
  states <- co$panel %>% count(wave,year,low_cr,state_dependency,name="n")
  write.csv(states,file.path(root,paste0("tables/",nm,"_wave_states.csv")),row.names=FALSE)
  missing <- co$baseline %>% summarise(across(c(age,female,education,bmi,smoking,hypertension,diabetes,heart,stroke,
     scr,cys,egfr_cr,egfr_cys,egfr_combined,blood_wt,severity),~sum(is.na(.x)))) %>%
    pivot_longer(everything(),names_to="variable",values_to="n_missing") %>% mutate(n_total=nrow(co$baseline),fraction=n_missing/n_total)
  write.csv(missing,file.path(root,paste0("tables/",nm,"_baseline_missingness.csv")),row.names=FALSE)
}
write.csv(bind_rows(audit),file.path(root,"audit/charls_linkage_and_flow.csv"),row.names=FALSE)
write.csv(bind_rows(raw_audit),file.path(root,"audit/charls_adl_coding_audit.csv"),row.names=FALSE)
distributions <- bind_rows(lapply(names(cohorts),function(nm) {
  d <- cohorts[[nm]]$baseline
  bind_rows(lapply(c("scr","cys","egfr_cr","egfr_cys","egfr_combined","bmi","blood_wt"),function(v) {
    x <- d[[v]][is.finite(d[[v]])]
    q <- if(length(x)) quantile(x,c(0,.01,.25,.5,.75,.99,1),names=FALSE) else rep(NA_real_,7)
    data.frame(cohort=nm,variable=v,n=length(x),min=q[1],p01=q[2],p25=q[3],median=q[4],p75=q[5],p99=q[6],max=q[7])
  }))
}))
write.csv(distributions,file.path(root,"tables/charls_marker_distributions.csv"),row.names=FALSE)
saveRDS(cohorts,file.path(root,"data_internal/charls_cohorts.rds"))
hashes <- data.frame(source=used,md5=unname(tools::md5sum(used)),bytes=file.info(used)$size)
write.csv(hashes,file.path(root,"audit/charls_used_source_hashes.csv"),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(root,"logs/R_session_build.txt"))
print(bind_rows(audit) %>% filter(grepl("baseline|cohort|panel",stage)))
message("Completed CHARLS panel. Only aggregate tables are exported; internal RDS has project-only row numbers.")

# Private linkage keys remain in the authorised raw workspace only.
private <- cohorts
private$baseline2011$baseline <- id_baselines[["1"]]
private$landmark2015$baseline <- id_baselines[["3"]]
dir.create(file.path(root,"data_internal/private_linkage"),recursive=TRUE,showWarnings=FALSE)
saveRDS(private,file.path(root,"data_internal/private_linkage/charls_baseline_keys.rds"))
