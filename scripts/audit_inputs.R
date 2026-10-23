suppressWarnings(Sys.setlocale("LC_CTYPE","English_United States.utf8"))
suppressPackageStartupMessages(library(jsonlite))
package <- normalizePath(commandArgs(TRUE)[1],winslash="/")
files <- list.files(package,pattern="[.]R$",recursive=TRUE,full.names=TRUE)
for(f in files) parse(f)
forbidden <- c("HHID","PN","ID","ID_w1","join_id","hhidpn","rahhidpn")
frames <- 0L; values <- 0L
walk <- function(x) {
  if(is.environment(x) || is.function(x)) stop("An executable environment was found in an input")
  if(is.data.frame(x)) {
    stopifnot(!any(names(x) %in% forbidden))
    frames <<- frames+1L
  }
  if(is.character(x)) {
    stopifnot(!any(grepl("[A-Za-z]:[/\\\\]",x)))
    values <<- values+length(x)
  }
  if(is.list(x)) for(y in x) walk(y)
}
for(name in c("clinical_panels.rds","renal_intervals.rds","selection_intervals.rds")) {
  file <- file.path(package,"local_data",name)
  if(file.exists(file)) walk(readRDS(file))
}
walk(readRDS(file.path(package,"data/bulk_expression_inputs.rds")))
report <- list(status="passed",R_scripts_parsed=length(files),data_frames_checked=frames,
  string_values_checked=values,raw_identifier_columns_found=FALSE,input_environments_found=FALSE,
  absolute_computer_paths_found=FALSE,scope="Clinical RDS observations and public bulk expression/design inputs")
write_json(report,file.path(package,"results/input_audit.json"),pretty=TRUE,auto_unbox=TRUE)
cat(toJSON(report,pretty=TRUE,auto_unbox=TRUE),"\n")
