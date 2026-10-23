args <- commandArgs(TRUE)
lib <- if(length(args)) args[1] else file.path(path.expand("~"),"R","review-reproduction-library")
dir.create(lib,recursive=TRUE,showWarnings=FALSE)
.libPaths(c(lib,.libPaths()))
options(repos=c(CRAN="https://cloud.r-project.org"))
packages <- c("survey","dplyr","tidyr","jsonlite","haven","readr","ggplot2","metafor","BiocManager")
missing <- packages[!vapply(packages,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing)) install.packages(missing,lib=lib)
bio <- c("limma","edgeR","statmod")
missing <- bio[!vapply(bio,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing)) BiocManager::install(missing,lib=lib,ask=FALSE,update=FALSE)
cat("R dependencies ready. Set REVIEW_R_LIB to: ",lib,"\n",sep="")
