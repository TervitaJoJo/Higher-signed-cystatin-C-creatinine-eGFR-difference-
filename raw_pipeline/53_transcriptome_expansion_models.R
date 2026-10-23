args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) %in% c(1L, 2L))
root <- normalizePath(args[1], winslash = "/")
out <- file.path(root, "transcriptome_expansion_20261006")
input <- file.path(out, "analysis_ready")

suppressPackageStartupMessages({
  library(limma)
  library(edgeR)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(jsonlite)
  library(data.table)
})
set.seed(20261006)
options(stringsAsFactors = FALSE)
for (folder in c("models", "qc", "tables", "figures")) {
  dir.create(file.path(out, folder), recursive = TRUE, showWarnings = FALSE)
}
sink(file.path(out, "model_run.log"), split = TRUE)
contracts <- fromJSON(file.path(input, "analysis_contracts.json"))
resume <- length(args) == 2L && args[2] == "--resume"
summaries <- if (resume) fromJSON(file.path(out, "model_summary.json"), simplifyVector = FALSE) else list()
annotations <- list()
key_cache <- new.env()
unique_value <- function(x) {
  x <- unique(x[!is.na(x) & nzchar(x)])
  if (length(x) == 1L) x else NA_character_
}
map_unique <- function(ids, type) {
  if (!exists(type, key_cache, inherits = FALSE)) {
    assign(type, keys(org.Hs.eg.db, keytype = type), key_cache)
  }
  valid <- intersect(unique(ids[!is.na(ids)]), get(type, key_cache))
  if (!length(valid)) return(rep(NA_character_, length(ids)))
  mapped <- suppressMessages(mapIds(org.Hs.eg.db, keys = valid, keytype = type,
                                   column = "SYMBOL", multiVals = unique_value))
  unname(mapped[match(ids, names(mapped))])
}
map_symbols <- function(ids) {
  valid <- keys(org.Hs.eg.db, keytype = "SYMBOL")
  mapped <- ifelse(ids %in% valid, ids, NA_character_)
  missing <- which(is.na(mapped) & !is.na(ids) & nzchar(ids))
  if (length(missing)) mapped[missing] <- map_unique(ids[missing], "ALIAS")
  mapped
}
map_list <- function(items, type) {
  values <- unique(unlist(items))
  mapped <- if (type == "SYMBOL") map_symbols(values) else map_unique(values, type)
  names(mapped) <- values
  vapply(items, function(x) unique_value(unname(mapped[x])), character(1))
}
split_identifiers <- function(ids) {
  lapply(strsplit(ifelse(is.na(ids), "", ids), "///", fixed = TRUE), function(x) unique(trimws(x[nzchar(x)])))
}
platform_annotation <- function(platform) {
  if (!is.null(annotations[[platform]])) return(annotations[[platform]])
  ann <- as.data.frame(fread(file.path(input, paste0(platform, "_annotation.csv")),
                             check.names = FALSE, colClasses = "character", encoding = "UTF-8"))
  ann$symbol <- map_list(split_identifiers(ann$source_entrez), "ENTREZID")
  from_symbol <- map_list(split_identifiers(ann$source_symbol), "SYMBOL")
  missing <- is.na(ann$symbol)
  ann$symbol[missing] <- from_symbol[missing]
  conflict <- !is.na(ann$symbol) & !is.na(from_symbol) & ann$symbol != from_symbol
  ann$symbol[conflict] <- NA_character_
  if (platform == "GPL16686") {
    refs <- lapply(ann$source_refseq, function(s) {
      if (is.na(s)) return(character())
      unique(regmatches(s, gregexpr("[NX][MR]_[0-9]+", s))[[1]])
    })
    ann$symbol <- map_list(refs, "REFSEQ")
  }
  stopifnot(!anyDuplicated(ann$feature_id))
  write.csv(ann, file.path(out, "qc", paste0(platform, "_mapped_annotation.csv")), row.names = FALSE)
  annotations[[platform]] <<- ann
  ann
}
load_data <- function(name) {
  cat("LOAD", name, "\n")
  contract <- contracts[contracts$name == name, ]
  stopifnot(nrow(contract) == 1L, !contract$provisional_only)
  tab <- read.delim(gzfile(file.path(input, paste0(name, "_expression.tsv.gz"))), check.names = FALSE)
  ids <- as.character(tab[[1]])
  x <- as.matrix(tab[-1])
  storage.mode(x) <- "double"
  rownames(x) <- ids
  m <- as.data.frame(fread(file.path(input, paste0(name, "_sample_sheet.csv")), check.names = FALSE, encoding = "UTF-8"))
  stopifnot(identical(colnames(x), m$matrix_column), !anyDuplicated(ids), !anyDuplicated(m$sample_accession), all(is.finite(x)))
  special <- grepl("^__", ids)
  if (any(special)) {
    write.csv(data.frame(feature_id = ids[special], x[special, , drop = FALSE], check.names = FALSE),
              file.path(out, "qc", paste0(name, "_excluded_htseq_rows.csv")), row.names = FALSE)
    x <- x[!special, , drop = FALSE]
    ids <- ids[!special]
  }
  a <- data.frame(feature_id = ids)
  if (contract$keytype == "ENSEMBL") {
    a$symbol <- map_unique(sub("\\..*$", "", ids), "ENSEMBL")
  } else if (contract$keytype == "SYMBOL") {
    a$symbol <- map_symbols(ids)
  } else if (contract$keytype == "SYMBOL_ENTREZ") {
    with_entrez <- grepl("\\|[0-9]+$", ids)
    a$symbol <- rep(NA_character_, length(ids))
    a$symbol[with_entrez] <- map_unique(sub("^.*\\|", "", ids[with_entrez]), "ENTREZID")
    a$symbol[!with_entrez] <- map_symbols(ids[!with_entrez])
  } else if (contract$keytype == "ENTREZ_SUFFIX") {
    ann <- platform_annotation(contract$platform)
    stopifnot(all(ids %in% ann$feature_id), all(grepl("^[0-9]+_at$", ids)))
    a$symbol <- map_unique(sub("_at$", "", ids), "ENTREZID")
  } else {
    ann <- platform_annotation(contract$platform)
    stopifnot(all(ids %in% ann$feature_id))
    a$symbol <- ann$symbol[match(ids, ann$feature_id)]
  }
  write.csv(a, file.path(out, "qc", paste0(name, "_feature_annotation.csv")), row.names = FALSE)
  for (field in intersect(c("age", "bmi"), names(m))) {
    m[[field]] <- suppressWarnings(as.numeric(m[[field]]))
  }
  list(x = x, m = m, a = a, contract = contract)
}
record_qc <- function(e, m, name, dge = NULL) {
  expression <- if (is.list(e)) e$E else e
  idx <- head(order(apply(expression, 1L, var), decreasing = TRUE), 2000L)
  pca <- prcomp(t(expression[idx, , drop = FALSE]), center = TRUE, scale. = FALSE)
  write.csv(data.frame(matrix_column = colnames(expression), PC1 = pca$x[, 1], PC2 = pca$x[, 2], m, check.names = FALSE),
            file.path(out, "qc", paste0(name, "_pca.csv")), row.names = FALSE)
  write.csv(data.frame(PC = seq_along(pca$sdev), variance_fraction = pca$sdev^2 / sum(pca$sdev^2)),
            file.path(out, "qc", paste0(name, "_pca_variance.csv")), row.names = FALSE)
  write.csv(cor(expression), file.path(out, "qc", paste0(name, "_sample_correlation.csv")))
  if (!is.null(dge)) write.csv(data.frame(matrix_column = colnames(expression), dge$samples),
                              file.path(out, "qc", paste0(name, "_library_sizes.csv")), row.names = FALSE)
}
fit_comparison <- function(d, reference, case, covariates, suffix, role, selection = rep(TRUE, nrow(d$m))) {
  chosen <- selection & d$m$group %in% c(reference, case)
  keep_sample <- complete.cases(d$m[, covariates, drop = FALSE])
  chosen <- chosen & keep_sample
  m <- d$m[chosen, , drop = FALSE]
  stopifnot(sum(m$group == reference) >= 3L, sum(m$group == case) >= 3L)
  if ("subject_id" %in% names(m)) stopifnot(!anyDuplicated(m$subject_id))
  m$group <- factor(m$group, levels = c(reference, case))
  for (field in covariates) {
    if (is.numeric(m[[field]])) m[[field]] <- m[[field]] - mean(m[[field]]) else m[[field]] <- factor(m[[field]])
  }
  design <- model.matrix(reformulate(c("group", covariates)), m)
  rownames(design) <- m$matrix_column
  name <- paste0(d$contract$name, "_", suffix)
  if (resume && name %in% names(summaries)) {
    stopifnot(file.exists(file.path(out, "tables", paste0(name, "_gene_de.csv"))),
              file.exists(file.path(out, "models", paste0(name, ".rds"))))
    cat("EXISTING_COMPLETE", name, "\n")
    return(invisible(NULL))
  }
  stopifnot(qr(design)$rank == ncol(design), nrow(design) > ncol(design) + 1L, all(is.finite(design)))
  write.csv(data.frame(matrix_column = rownames(design), design, check.names = FALSE),
            file.path(out, "qc", paste0(name, "_design.csv")), row.names = FALSE)
  x <- d$x[, chosen, drop = FALSE]
  if (d$contract$assay == "counts") {
    stopifnot(all(x >= 0))
    y <- DGEList(x)
    retained <- filterByExpr(y, design = design)
    stopifnot(sum(retained) >= 500L)
    y <- calcNormFactors(y[retained, , keep.lib.sizes = FALSE], method = "TMM")
    pdf(file.path(out, "figures", paste0(name, "_voom.pdf")), width = 6, height = 5)
    e <- voom(y, design, plot = TRUE)
    dev.off()
    trend <- FALSE
    assay <- "TMM_limma_voom"
    record_qc(e, m, name, y)
  } else {
    if (d$contract$assay == "array_MAS5_linear") {
      x <- log2(x + 1)
      assay <- "MAS5_log2_intensity_plus1_limma_trend"
    } else assay <- "deposited_log_scale_array_limma_trend"
    retained <- apply(x, 1L, var) > 0
    e <- x[retained, , drop = FALSE]
    trend <- TRUE
    record_qc(e, m, name)
  }
  write.csv(data.frame(feature_id = rownames(x), retained = retained),
            file.path(out, "qc", paste0(name, "_expression_filter.csv")), row.names = FALSE)
  fit <- eBayes(lmFit(e, design), trend = trend, robust = TRUE)
  coefficient <- paste0("group", case)
  index <- match(coefficient, colnames(design))
  stopifnot(!is.na(index))
  tt <- topTable(fit, coef = index, number = Inf, sort.by = "none", confint = TRUE)
  stopifnot(identical(rownames(tt), d$a$feature_id[retained]))
  tt <- cbind(d$a[retained, ], tt)
  tt$SE <- fit$stdev.unscaled[, index] * sqrt(fit$s2.post)
  tt$df_total <- fit$df.total
  tt$contrast <- name
  tt$direction <- paste(case, "minus", reference)
  valid <- which(!is.na(tt$symbol) & nzchar(tt$symbol))
  representatives <- valid[order(-tt$AveExpr[valid], tt$feature_id[valid])]
  representatives <- representatives[!duplicated(tt$symbol[representatives])]
  genes <- tt[representatives, , drop = FALSE]
  genes$gene_FDR <- p.adjust(genes$P.Value, "BH")
  genes <- genes[order(genes$gene_FDR, genes$P.Value), ]
  write.csv(tt, file.path(out, "tables", paste0(name, "_feature_de.csv")), row.names = FALSE)
  write.csv(genes, file.path(out, "tables", paste0(name, "_gene_de.csv")), row.names = FALSE)
  saveRDS(list(e = e, m = m, a = d$a[retained, ], design = design, fit = fit, coefficient = coefficient),
          file.path(out, "models", paste0(name, ".rds")))
  summaries[[name]] <<- list(contrast = name, accession = d$contract$accession, role = role,
      direction = paste(case, "minus", reference), assay = assay, covariates = covariates,
      n = nrow(m), group_counts = as.list(table(m$group)), excluded_missing_covariates = sum(selection & d$m$group %in% c(reference, case) & !keep_sample),
      features_input = nrow(d$x), features_tested = nrow(tt), genes_tested = nrow(genes),
      feature_fdr05 = sum(tt$adj.P.Val < .05), gene_fdr05 = sum(genes$gene_FDR < .05),
      strict_gene_and_feature_fdr05 = sum(genes$gene_FDR < .05 & genes$adj.P.Val < .05),
      rank = qr(design)$rank, residual_df = unique(fit$df.residual), condition_number = kappa(design))
  write_json(summaries, file.path(out, "model_summary.json"), pretty = TRUE, auto_unbox = TRUE, na = "null")
  cat("COMPLETE", name, "n", nrow(m), "features", nrow(tt), "gene_FDR05", sum(genes$gene_FDR < .05), "\n")
}

for (name in c("GSE111006", "GSE111010", "GSE111016")) {
  d <- load_data(name)
  fit_comparison(d, "Healthy", "Sarcopenia", "age", "sarcopenia", "muscle_function_primary")
}
d <- load_data("GSE117525")
stopifnot(sum(d$m$group %in% c("Healthy", "Frail")) == 134L)
fit_comparison(d, "Healthy", "Frail", c("age", "sex", "bmi"), "frailty_baseline", "muscle_function_primary")
d$m$plate_proxy <- sub("_.*$", "", d$m$title)
stopifnot(all(d$m$plate_proxy %in% c("G084", "G085", "G086", "G087")))
fit_comparison(d, "Healthy", "Frail", c("age", "sex", "bmi", "plate_proxy"),
               "frailty_title_prefix", "muscle_function_title_prefix_sensitivity")
d <- load_data("GSE174106")
fit_comparison(d, "Stablehigh", "ALMloss", c("age", "sex"), "alm_loss", "muscle_mass_secondary")
d <- load_data("GSE129534")
fit_comparison(d, "Nonfrail", "Frail", c("age", "sex", "race"), "frailty", "middle_age_pbmc_secondary")
d <- load_data("GSE140358")
fit_comparison(d, "Robust", "Frail", "age", "frailty", "older_wholeblood_primary")
d <- load_data("GSE65218")
fit_comparison(d, "Nonfrail", "Frail", c("sex", "bmi"), "frailty", "older_pbmc_primary")
fit_comparison(d, "Nonfrail", "Prefrail", c("sex", "bmi"), "prefrailty", "older_pbmc_secondary")
d <- load_data("GSE37171")
stopifnot(sum(d$m$study_set == "discovery") == 83L)
fit_comparison(d, "Healthy", "Uremia", c("age", "sex"), "uremia_discovery", "renal_wholeblood_primary", d$m$study_set == "discovery")
fit_comparison(d, "Healthy", "Uremia", c("age", "sex", "race"), "uremia_race", "renal_wholeblood_sensitivity", d$m$study_set == "discovery")
d <- load_data("GSE15072_GPL96")
fit_comparison(d, "Healthy", "CKD", character(), "ckd", "renal_pbmc_primary")
fit_comparison(d, "Healthy", "HD", character(), "hemodialysis", "renal_pbmc_secondary")
d <- load_data("GSE142153")
fit_comparison(d, "Healthy", "DN", character(), "dn", "renal_pbmc_diabetes_support")
fit_comparison(d, "Healthy", "ESRD", character(), "esrd", "renal_pbmc_diabetes_support")
capture.output(sessionInfo(), file = file.path(out, "analysis_sessionInfo.txt"))
cat("ALL_EXPANSION_MODELS_COMPLETE\n")
sink()
