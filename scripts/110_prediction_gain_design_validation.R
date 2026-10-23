Sys.setlocale("LC_CTYPE", "English_United States.utf8")
suppressPackageStartupMessages({library(dplyr); library(jsonlite)})
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args) >= 1L)
root <- normalizePath(args[1], winslash="/", mustWork=TRUE)
B <- if (length(args) > 1L) as.integer(args[2]) else 300L
out <- file.path(root, "prediction_gain_20261008", "design_validation")
for (folder in c("tables", "audit", "reports", "data_internal")) dir.create(file.path(out,folder),recursive=TRUE,showWarnings=FALSE)
cache <- readRDS(file.path(root,"prediction_gain_20261008","data_internal_prediction_cache.rds"))

auc <- function(y,p,w) {
  ord <- order(p); y <- y[ord]; w <- w[ord]; p <- p[ord]
  grp <- cumsum(c(TRUE,diff(p)!=0))
  wp <- as.numeric(rowsum(w*(y==1),grp,reorder=FALSE))
  wn <- as.numeric(rowsum(w*(y==0),grp,reorder=FALSE))
  sum(wp*(cumsum(wn)-wn/2))/(sum(wp)*sum(wn))
}
wm <- function(x,w) sum(x*w)/sum(w)
metrics <- function(y,p0,p1,w) {
  ok <- is.finite(p0)&is.finite(p1)&is.finite(w)&w>0
  y<-y[ok]; p0<-p0[ok]; p1<-p1[ok]; w<-w[ok]
  if (sum(y==1)<10 || sum(y==0)<10) return(rep(NA_real_,13))
  ev<-y==1; no<-y==0
  nri_e<-wm(sign(p1[ev]-p0[ev]),w[ev])
  nri_n<-wm(sign(p0[no]-p1[no]),w[no])
  a0<-auc(y,p0,w); a1<-auc(y,p1,w)
  b0<-wm((p0-y)^2,w); b1<-wm((p1-y)^2,w)
  cal <- function(p) {
    lp<-qlogis(pmin(pmax(p,1e-6),1-1e-6))
    intercept<-coef(glm(y~1+offset(lp),family=quasibinomial(),weights=w))[1]
    slope<-coef(glm(y~lp,family=quasibinomial(),weights=w))[2]
    c(intercept,slope)
  }
  c(a0,a1,a1-a0,b0,b1,b1-b0,nri_e+nri_n,nri_e,nri_n,
    wm(p1[ev]-p0[ev],w[ev])-wm(p1[no]-p0[no],w[no]),cal(p0)[1],cal(p1)[1],cal(p1)[2])
}
metric_names<-c("base_auc","plus_auc","delta_auc","base_brier","plus_brier","delta_brier",
  "continuous_nri","nri_event","nri_nonevent","idi","base_calibration_intercept","plus_calibration_intercept","plus_calibration_slope")

crossfit <- function(d,terms,kidney,w,fold) {
  p0<-p1<-rep(NA_real_,nrow(d))
  f0<-reformulate(c(terms,kidney),response="recovery")
  f1<-reformulate(c(terms,kidney,"delta10"),response="recovery")
  for (f in 1:5) {
    tr<-fold!=f&w>0; te<-fold==f&w>0
    if (!any(te) || sum(d$recovery[tr]==1)<20 || sum(d$recovery[tr]==0)<20) next
    train<-d[tr,,drop=FALSE]; train$fit_weight<-w[tr]/mean(w[tr])
    fits<-tryCatch({
      m0<-glm(f0,data=train,weights=fit_weight,family=quasibinomial())
      m1<-glm(f1,data=train,weights=fit_weight,family=quasibinomial())
      stopifnot(m0$converged,m1$converged)
      list(predict(m0,newdata=d[te,,drop=FALSE],type="response"),predict(m1,newdata=d[te,,drop=FALSE],type="response"))
    },error=function(e)NULL)
    if (!is.null(fits)) {p0[te]<-fits[[1]];p1[te]<-fits[[2]]}
  }
  if (any(!is.finite(p0[w>0])) || any(!is.finite(p1[w>0]))) return(NULL)
  list(p0=p0,p1=p1,metrics=metrics(d$recovery,p0,p1,w))
}

resample_psu <- function(psus) {
  mult<-rep(0,nrow(psus))
  for (ii in split(seq_len(nrow(psus)),psus$stratum)) {
    nh<-length(ii)
    if (nh==1L) mult[ii]<-1 else mult[ii]<-tabulate(sample.int(nh,nh-1L,replace=TRUE),nbins=nh)*nh/(nh-1L)
  }
  mult
}

rows<-list(); saved<-list(); calibration<-list(); bootall<-list()
for (cohort in names(cache)) {
  d<-as.data.frame(cache[[cohort]])
  d$cluster_key<-interaction(d$stratum,d$psu,drop=TRUE)
  psus<-distinct(d,cluster_key,stratum)
  pi<-match(d$cluster_key,psus$cluster_key)
  terms<-c("age_c","female","education","interval","start_severity","bmi","smoking","hypertension","diabetes","heart","stroke")
  if (cohort=="HRS2016VBS") terms<-c(terms,"draw_lag_months")
  terms<-terms[terms%in%names(d)]
  terms<-terms[vapply(d[terms],function(x)length(unique(x))>1L,logical(1))]
  set.seed(20261008+nrow(d))
  fold_psu<-sample(rep(1:5,length.out=nrow(psus)))
  fold<-fold_psu[pi]
  weights<-d$analysis_weight
  for (kidney in c("mean_egfr","egfr_cr","egfr_cys")) {
    point<-crossfit(d,terms,kidney,weights,fold)
    stopifnot(!is.null(point))
    vals<-matrix(NA_real_,B,length(metric_names),dimnames=list(NULL,metric_names))
    set.seed(20261009+nrow(d))
    for (b in seq_len(B)) {
      mult<-resample_psu(psus)
      wb<-weights*mult[pi]
      fold_b<-sample(rep(1:5,length.out=nrow(psus)))[pi]
      fit<-crossfit(d,terms,kidney,wb,fold_b)
      if (!is.null(fit)) vals[b,]<-fit$metrics
      if (b%%50L==0L) cat(cohort,kidney,"bootstrap",b,"/",B,"\n")
    }
    for (j in seq_along(metric_names)) {
      ok<-is.finite(vals[,j]); ci<-quantile(vals[ok,j],c(.025,.975),na.rm=TRUE)
      rows[[length(rows)+1L]]<-data.frame(cohort=cohort,kidney_base=kidney,metric=metric_names[j],
        estimate=point$metrics[j],ci_low=ci[1],ci_high=ci[2],bootstrap_valid=sum(ok),bootstrap_requested=B,
        n_people=nrow(d),n_psu=nrow(psus),events=sum(d$recovery==1))
    }
    saved[[paste(cohort,kidney,sep="_")]]<-data.frame(pid=d$pid,recovery=d$recovery,weight=weights,
      cluster=d$cluster_key,fold=fold,pred_base=point$p0,pred_plus=point$p1)
    bootall[[paste(cohort,kidney,sep="_")]]<-vals
    for (model in c("base","plus")) {
      p<-if(model=="base")point$p0 else point$p1
      bin<-pmin(10L,ceiling(rank(p,ties.method="average")/length(p)*10))
      for (ii in split(seq_along(p),bin)) calibration[[length(calibration)+1L]]<-data.frame(
        cohort=cohort,kidney_base=kidney,model=model,bin=bin[ii[1]],n=length(ii),
        predicted=wm(p[ii],weights[ii]),observed=wm(d$recovery[ii],weights[ii]))
    }
  }
}
res<-bind_rows(rows)
write.csv(res,file.path(out,"tables","prediction_design_validated_metrics.csv"),row.names=FALSE)
write.csv(bind_rows(calibration),file.path(out,"tables","calibration_deciles.csv"),row.names=FALSE)
saveRDS(saved,file.path(out,"data_internal","out_of_fold_predictions.rds"))
saveRDS(bootall,file.path(out,"data_internal","refit_bootstrap_metrics.rds"))
write_json(list(outcome="Recovery at the first eligible next visit; death included as non-recovery",
  sampling="Complete-case participants from the existing common clinical sample",
  folds="Five folds at survey PSU level, keeping every participant from a PSU in one fold",
  uncertainty="Rao-Wu rescaled PSU bootstrap within strata, refitting and cross-validating both models in every replicate",
  singleton_strata="One PSU retained with multiplier 1; uncertainty conditional on singleton PSUs",
  weight_model="Weighted quasibinomial GLM; fitted means equivalent to single-level svyglm",
  estimand="Exploratory incremental next-visit prediction; interval varies across cohorts and participants",
  base_models=c("clinical covariates + mean eGFR","clinical covariates + eGFRcr","clinical covariates + eGFRcys"),
  plus="Each base + delta10; do not add delta10 to a model already containing both kidney markers",
  continuous_nri="Directional event and nonevent components; not percentage correctly reclassified",
  limitations=c("Internal validation, no independent external predictive validation","No fixed-time clinical deployment model","Outcome-observed complete-case selection remains","Few PSUs per HRS stratum limit effective inference"),
  bootstrap_replicates=B),file.path(out,"audit","prediction_design_validation_specification.json"),pretty=TRUE,auto_unbox=TRUE)
capture.output(sessionInfo(),file=file.path(out,"audit","R_sessionInfo.txt"))
cat("PREDICTION_DESIGN_VALIDATION_COMPLETE\n")
print(res%>%filter(metric%in%c("delta_auc","continuous_nri","delta_brier","plus_calibration_slope")))
