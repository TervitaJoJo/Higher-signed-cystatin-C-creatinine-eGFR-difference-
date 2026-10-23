egfr_cr <- function(scr,age,female) {
  k <- ifelse(female==1,.7,.9);a <- ifelse(female==1,-.241,-.302)
  out <- 142*pmin(scr/k,1)^a*pmax(scr/k,1)^(-1.200)*.9938^age*ifelse(female==1,1.012,1)
  ifelse(is.finite(out)&scr>0&age>0&!is.na(female),out,NA_real_)
}
egfr_cys <- function(cys,age,female) {
  out <- 133*pmin(cys/.8,1)^(-.499)*pmax(cys/.8,1)^(-1.328)*.996^age*ifelse(female==1,.932,1)
  ifelse(is.finite(out)&cys>0&age>0&!is.na(female),out,NA_real_)
}
egfr_combined <- function(scr,cys,age,female) {
  k <- ifelse(female==1,.7,.9);a <- ifelse(female==1,-.219,-.144)
  out <- 135*pmin(scr/k,1)^a*pmax(scr/k,1)^(-.544)*pmin(cys/.8,1)^(-.323)*
    pmax(cys/.8,1)^(-.778)*.9961^age*ifelse(female==1,.963,1)
  ifelse(is.finite(out)&scr>0&cys>0&age>0&!is.na(female),out,NA_real_)
}
stopifnot(abs(egfr_cr(.9,60,0)-97.768)<.01,abs(egfr_cys(.8,60,0)-104.57)<.05,
  abs(egfr_cr(.7,60,1)-142*.9938^60*1.012)<1e-10,
  abs(egfr_combined(.9,.8,60,0)-135*.9961^60)<1e-10,
  is.na(egfr_cr(-1,60,0)),is.na(egfr_combined(1,NA,60,1)),is.na(egfr_cys(1,60,NA)))
