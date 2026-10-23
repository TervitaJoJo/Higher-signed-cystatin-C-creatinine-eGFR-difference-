args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
root <- normalizePath(args[1], winslash = "/")
out <- file.path(root, "transcriptome_analysis")

suppressPackageStartupMessages({
  library(limma)
  library(edgeR)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(jsonlite)
})
set.seed(20261006)
options(stringsAsFactors = FALSE)
input <- file.path(root, "transcriptome_feasibility", "analysis_ready")
resources <- file.path(out, "resources")
for (folder in c("models", "qc", "tables", "figures")) {
  dir.create(file.path(out, folder), recursive = TRUE, showWarnings = FALSE)
}
log_file <- file.path(out, "model_run.log")
sink(log_file, split = TRUE)
summaries <- list()
input_audits <- list()
libraries <- c("MSigDB_Hallmark_2020", "Reactome_2022")
sets <- unlist(lapply(libraries, function(lib) {
  lines <- strsplit(readLines(file.path(resources, paste0(lib, ".gmt")), warn = FALSE), "\t")
  paths <- lapply(lines, function(x) unique(x[-c(1, 2)]))
  names(paths) <- paste(lib, vapply(lines, `[`, character(1), 1), sep = "::")
  paths
}), recursive = FALSE)
stopifnot(!anyDuplicated(names(sets)))

unique_map <- function(x) {
  x <- unique(x[!is.na(x) & nzchar(x)])
  if (length(x) == 1L) x else NA_character_
}
map_unique <- function(ids, keytype, column) {
  valid <- intersect(unique(ids[!is.na(ids)]), keys(org.Hs.eg.db, keytype = keytype))
  if (!length(valid)) return(rep(NA_character_, length(ids)))
  mapped <- suppressMessages(mapIds(org.Hs.eg.db, keys = valid, column = column,
                                    keytype = keytype, multiVals = unique_map))
  unname(mapped[match(ids, names(mapped))])
}

load_data <- function(accession, probe = FALSE) {
  cat("\nLOAD", accession, "\n")
  tab <- read.delim(gzfile(file.path(input, paste0(accession, "_expression.tsv.gz"))),
                    check.names = FALSE)
  ids <- as.character(tab[[1]])
  x <- as.matrix(tab[-1])
  storage.mode(x) <- "double"
  rownames(x) <- ids
  m <- read.csv(file.path(input, paste0(accession, "_sample_sheet.csv")),
                check.names = FALSE, fileEncoding = "UTF-8-BOM")
  stopifnot("sample_accession" %in% names(m), length(m$sample_accession) == ncol(x))
  stopifnot(identical(colnames(x), m$matrix_column), !anyDuplicated(ids),
            !anyDuplicated(m$sample_accession), all(is.finite(x)))
  raw_features <- nrow(x)
  special <- grepl("^__", ids)
  if (any(special)) {
    write.csv(data.frame(feature_id = ids[special], x[special, , drop = FALSE], check.names = FALSE),
              file.path(out, "qc", paste0(accession, "_excluded_htseq_summary_rows.csv")), row.names = FALSE)
    x <- x[!special, , drop = FALSE]
    ids <- ids[!special]
  }
  annotation <- data.frame(feature_id = ids)
  if (probe) {
    ann <- read.csv(file.path(resources, "GPL10558_probe_annotation.csv"), check.names = FALSE)
    idx <- match(ids, ann$feature_id)
    stopifnot(!anyNA(idx))
    entrez <- as.character(ann$entrez_id[idx])
    entrez[is.na(entrez) | !grepl("^[0-9]+$", entrez)] <- NA_character_
    annotation$entrez_id <- entrez
    annotation$source_symbol <- ann$source_symbol[idx]
    annotation$symbol <- map_unique(entrez, "ENTREZID", "SYMBOL")
    annotation$ensembl <- map_unique(entrez, "ENTREZID", "ENSEMBL")
  } else {
    annotation$ensembl <- sub("[._].*$", "", ids)
    stopifnot(all(grepl("^ENSG[0-9]+$", annotation$ensembl)))
    annotation$symbol <- map_unique(annotation$ensembl, "ENSEMBL", "SYMBOL")
  }
  write.csv(annotation, file.path(out, "qc", paste0(accession, "_feature_annotation.csv")), row.names = FALSE)
  input_audits[[accession]] <<- list(raw_features = raw_features, features = nrow(x),
      excluded_htseq_summary_rows = sum(special), samples = ncol(x),
      unique_symbols = length(unique(na.omit(annotation$symbol))),
      annotated_features = sum(!is.na(annotation$symbol)),
      fractional_entry_fraction = mean(abs(x - round(x)) > 1e-7),
      negative_entries = sum(x < 0))
  list(x = x, m = m, a = annotation, accession = accession)
}

assert_design <- function(design, name) {
  stopifnot(all(is.finite(design)), qr(design)$rank == ncol(design),
            nrow(design) > ncol(design) + 1L)
  write.csv(data.frame(matrix_column = rownames(design), design, check.names = FALSE),
            file.path(out, "qc", paste0(name, "_design.csv")), row.names = FALSE)
}

record_qc <- function(e, metadata, name, dge = NULL) {
  matrix <- if (is.list(e)) e$E else e
  variance <- apply(matrix, 1, var)
  idx <- head(order(variance, decreasing = TRUE), 2000L)
  pca <- prcomp(t(matrix[idx, , drop = FALSE]), center = TRUE, scale. = FALSE)
  fraction <- pca$sdev^2 / sum(pca$sdev^2)
  pcs <- data.frame(matrix_column = colnames(matrix), PC1 = pca$x[, 1], PC2 = pca$x[, 2],
                    metadata[, setdiff(names(metadata), "matrix_column"), drop = FALSE], check.names = FALSE)
  write.csv(pcs, file.path(out, "qc", paste0(name, "_pca.csv")), row.names = FALSE)
  write.csv(data.frame(PC = seq_along(fraction), variance_fraction = fraction),
            file.path(out, "qc", paste0(name, "_pca_variance.csv")), row.names = FALSE)
  write.csv(cor(matrix), file.path(out, "qc", paste0(name, "_sample_correlation.csv")))
  library <- data.frame(matrix_column = colnames(matrix), median_expression = apply(matrix, 2, median))
  if (!is.null(dge)) library <- cbind(library, dge$samples[, c("lib.size", "norm.factors")])
  write.csv(library, file.path(out, "qc", paste0(name, "_libraries.csv")), row.names = FALSE)
  saveRDS(list(expression = matrix, metadata = metadata, feature_ids = rownames(matrix)),
          file.path(out, "qc", paste0(name, "_normalized_expression.rds")))
}

count_input <- function(data, select, design, name, block = NULL) {
  x <- data$x[, select, drop = FALSE]
  stopifnot(all(x >= 0))
  d <- DGEList(counts = x)
  keep <- filterByExpr(d, design = design)
  stopifnot(sum(keep) > 1000L)
  d <- calcNormFactors(d[keep, , keep.lib.sizes = FALSE], method = "TMM")
  pdf(file.path(out, "figures", paste0(name, "_voom.pdf")), width = 6, height = 5)
  v <- voom(d, design, plot = TRUE)
  dev.off()
  correlation <- NULL
  if (!is.null(block)) {
    first <- duplicateCorrelation(v, design, block = block)$consensus.correlation
    v <- voom(d, design, block = block, correlation = first, plot = FALSE)
    correlation <- duplicateCorrelation(v, design, block = block)$consensus.correlation
    stopifnot(is.finite(correlation))
    v <- voom(d, design, block = block, correlation = correlation, plot = FALSE)
  }
  record_qc(v, data$m[select, , drop = FALSE], name, d)
  write.csv(data.frame(feature_id = rownames(x), retained = keep),
            file.path(out, "qc", paste0(name, "_expression_filter.csv")), row.names = FALSE)
  list(e = v, a = data$a[keep, , drop = FALSE], design = design,
       m = data$m[select, , drop = FALSE], correlation = correlation, block = block,
       assay = "TMM_limma_voom", trend = FALSE, input_features = nrow(x))
}

enrich <- function(model, contrast, table, name) {
  # One expression-based representative per gene prevents probe multiplicity bias.
  valid <- which(!is.na(table$symbol) & nzchar(table$symbol))
  idx <- valid[order(-table$AveExpr[valid], table$feature_id[valid])]
  idx <- idx[!duplicated(table$symbol[idx])]
  symbol <- table$symbol[idx]
  indices <- lapply(sets, function(s) which(symbol %in% s))
  indices <- indices[lengths(indices) >= 15 & lengths(indices) <= 500]
  stopifnot(length(indices) > 40L)
  if (is.null(model$block)) {
    expression <- if (is.list(model$e)) model$e$E[idx, , drop = FALSE] else model$e[idx, , drop = FALSE]
    weights <- if (is.list(model$e)) model$e$weights[idx, , drop = FALSE] else NULL
    result <- camera(expression, indices, design = model$design, contrast = contrast,
                     weights = weights, inter.gene.cor = 0.01, trend.var = model$trend, sort = FALSE)
    method <- "CAMERA_inter_gene_cor_0.01"
  } else {
    result <- cameraPR(table$t[idx], indices, inter.gene.cor = 0.01, sort = FALSE)
    method <- "CAMERA_PR_block_adjusted_t_inter_gene_cor_0.01"
  }
  result$pathway_id <- rownames(result)
  result$library <- sub("::.*$", "", result$pathway_id)
  result$term <- sub("^[^:]+::", "", result$pathway_id)
  result$FDR_cross_library <- p.adjust(result$PValue, "BH")
  result$FDR_within_library <- ave(result$PValue, result$library, FUN = function(p) p.adjust(p, "BH"))
  result$method <- method
  result$contrast <- name
  result <- result[order(result$FDR_cross_library, result$PValue), ]
  write.csv(result, file.path(out, "tables", paste0(name, "_camera.csv")), row.names = FALSE)
  ranked <- table[idx, c("feature_id", "ensembl", "symbol", "AveExpr", "logFC", "t", "P.Value", "adj.P.Val")]
  write.csv(ranked[order(-ranked$t), ], file.path(out, "tables", paste0(name, "_ranked_symbols.csv")), row.names = FALSE)
  write.csv(data.frame(pathway_id = names(indices), observed_genes = lengths(indices),
                      source_genes = lengths(sets[names(indices)])),
            file.path(out, "qc", paste0(name, "_pathway_coverage.csv")), row.names = FALSE)
  sum(result$FDR_cross_library < 0.05)
}

fit_model <- function(model, contrasts, name, labels, role, save_model = TRUE) {
  assert_design(model$design, name)
  if (is.null(model$block)) fit <- lmFit(model$e, model$design) else {
    fit <- lmFit(model$e, model$design, block = model$block, correlation = model$correlation)
  }
  fitted <- eBayes(contrasts.fit(fit, contrasts), trend = model$trend, robust = TRUE)
  if (save_model) saveRDS(list(model = model, fit = fitted, contrasts = contrasts),
                         file.path(out, "models", paste0(name, ".rds")))
  for (i in seq_len(ncol(contrasts))) {
    contrast_name <- colnames(contrasts)[i]
    tt <- topTable(fitted, coef = i, number = Inf, sort.by = "none", confint = TRUE)
    stopifnot(identical(rownames(tt), model$a$feature_id))
    tt <- cbind(model$a, tt)
    tt$contrast <- contrast_name
    tt$direction_definition <- labels[i]
    write.csv(tt, file.path(out, "tables", paste0(contrast_name, "_de.csv")), row.names = FALSE)
    enriched <- enrich(model, contrasts[, i], tt, contrast_name)
    summaries[[contrast_name]] <<- list(contrast = contrast_name, model = name, role = role,
        direction = labels[i], n = nrow(model$m), input_features = model$input_features,
        tested_features = nrow(tt), mapped_symbols = length(unique(na.omit(tt$symbol))),
        residual_df = unique(fitted$df.residual), design_rank = qr(model$design)$rank,
        condition_number = kappa(model$design), assay = model$assay,
        group_counts = if ("group" %in% names(model$m)) as.list(table(model$m$group)) else NULL,
        de_fdr05 = sum(tt$adj.P.Val < 0.05, na.rm = TRUE),
        de_fdr05_abs_logfc05 = sum(tt$adj.P.Val < 0.05 & abs(tt$logFC) >= 0.5, na.rm = TRUE),
        camera_cross_library_fdr05 = enriched, block_correlation = model$correlation)
    cat("COMPLETE", contrast_name, "n=", nrow(model$m), "tested=", nrow(tt),
        "DE_FDR05=", sum(tt$adj.P.Val < 0.05), "CAMERA_FDR05=", enriched, "\n")
    write_json(summaries, file.path(out, "model_summary.json"), pretty = TRUE, auto_unbox = TRUE, na = "null")
  }
  invisible(model)
}

single_contrast <- function(design, column, name) {
  contrast <- matrix(0, ncol(design), 1L, dimnames = list(colnames(design), name))
  contrast[column, 1] <- 1
  contrast
}

# Older-only functional comparison; no young samples or derived logit groups.
d <- load_data("GSE144304")
for (field in c("age", "bmi", "sppb", "meangrip", "walktime")) {
  d$m[[field]] <- suppressWarnings(as.numeric(d$m[[field]]))
}
sel <- d$m$group %in% c("Fit", "Frail")
m <- d$m[sel, ]
stopifnot(nrow(m) == 54L, all(complete.cases(m[, c("age", "gender", "bmi", "sppb", "meangrip", "walktime")])))
m$group <- factor(m$group, levels = c("Fit", "Frail"))
m$sex <- factor(m$gender)
m$age <- m$age - mean(m$age)
m$bmi <- m$bmi - mean(m$bmi)
design <- model.matrix(~ group + age + sex + bmi, m)
rownames(design) <- m$matrix_column
model <- count_input(d, sel, design, "GSE144304_older")
fit_model(model, single_contrast(design, "groupFrail", "GSE144304_frail_vs_fit"),
          "GSE144304_older", "Frail minus Fit; age/sex/BMI adjusted", "primary")
for (phenotype in c("sppb", "meangrip", "walktime")) {
  secondary <- m
  secondary$worse_function <- as.numeric(scale(secondary[[phenotype]])) * if (phenotype == "walktime") 1 else -1
  design <- model.matrix(~ worse_function + age + sex + bmi, secondary)
  rownames(design) <- secondary$matrix_column
  name <- paste0("GSE144304_worse_", phenotype)
  model <- count_input(d, sel, design, name)
  fit_model(model, single_contrast(design, "worse_function", name), name,
            paste0("Per SD worse ", phenotype, "; age/sex/BMI adjusted"), "secondary")
}

# RSEM expected counts stay numeric; independent donors form the primary model.
d <- load_data("GSE157712")
sel <- d$m$group %in% c("Control", "CKD")
m <- d$m[sel, ]
stopifnot(nrow(m) == 14L, !anyDuplicated(m$`subject id`))
m$group <- factor(m$group, levels = c("Control", "CKD"))
m$sex <- factor(m$sex)
m$age <- m$age - mean(m$age)
design <- model.matrix(~ group + age + sex, m)
rownames(design) <- m$matrix_column
model <- count_input(d, sel, design, "GSE157712_independent")
fit_model(model, single_contrast(design, "groupCKD", "GSE157712_ckd_vs_control"),
          "GSE157712_independent", "CKD minus Control; age/sex adjusted", "primary")
m$race <- factor(m$race)
design <- model.matrix(~ group + age + sex + race, m)
rownames(design) <- m$matrix_column
model <- count_input(d, sel, design, "GSE157712_race_sensitivity")
fit_model(model, single_contrast(design, "groupCKD", "GSE157712_ckd_race_adjusted"),
          "GSE157712_race_sensitivity", "CKD minus Control; age/sex/race adjusted", "sensitivity")
m <- d$m
m$group <- factor(m$group, levels = c("Control", "CKD", "Dialysis"))
m$sex <- factor(m$sex)
m$age <- m$age - mean(m$age)
design <- model.matrix(~ group + age + sex, m)
rownames(design) <- m$matrix_column
model <- count_input(d, rep(TRUE, nrow(m)), design, "GSE157712_subject_blocked", block = m$`subject id`)
contrasts <- cbind(single_contrast(design, "groupCKD", "GSE157712_blocked_ckd_vs_control"),
                   single_contrast(design, "groupDialysis", "GSE157712_blocked_dialysis_vs_control"))
fit_model(model, contrasts, "GSE157712_subject_blocked",
          c("CKD minus Control; subject blocked; age/sex adjusted",
            "Dialysis minus Control; subject blocked; age/sex adjusted"), "sensitivity")

# TPM supports only a processed-expression sensitivity analysis.
d <- load_data("GSE287118")
m <- d$m
m$group <- factor(m$`disease state`, levels = c("Control", "CKD", "PAD", "PAD+CKD"))
sel <- m$group %in% c("Control", "CKD")
design <- model.matrix(~ 0 + group, droplevels(m[sel, ]))
colnames(design) <- c("Control", "CKD")
rownames(design) <- m$matrix_column[sel]
keep <- rowSums(d$x[, sel, drop = FALSE] >= 1) >= 16L
e <- log2(d$x[keep, sel, drop = FALSE] + 0.5)
record_qc(e, m[sel, ], "GSE287118_no_pad")
model <- list(e = e, a = d$a[keep, ], design = design, m = m[sel, ],
              block = NULL, correlation = NULL, assay = "log2_TPM_plus_0.5_limma_trend",
              trend = TRUE, input_features = nrow(d$x))
fit_model(model, makeContrasts(GSE287118_ckd_no_pad = CKD - Control, levels = design),
          "GSE287118_no_pad", "CKD minus Control; no PAD; no demographic adjustment possible", "exploratory_support")
design <- model.matrix(~ 0 + group, m)
colnames(design) <- c("Control", "CKD", "PAD", "PADCKD")
rownames(design) <- m$matrix_column
keep <- rowSums(d$x >= 1) >= 14L
e <- log2(d$x[keep, ] + 0.5)
record_qc(e, m, "GSE287118_four_groups")
model <- list(e = e, a = d$a[keep, ], design = design, m = m,
              block = NULL, correlation = NULL, assay = "log2_TPM_plus_0.5_limma_trend",
              trend = TRUE, input_features = nrow(d$x))
contrasts <- makeContrasts(GSE287118_ckd_no_pad_four_groups = CKD - Control,
                         GSE287118_ckd_with_pad = PADCKD - PAD,
                         GSE287118_ckd_pad_interaction = (PADCKD - PAD) - (CKD - Control), levels = design)
fit_model(model, contrasts, "GSE287118_four_groups",
          c("CKD minus Control in four-group model", "PAD+CKD minus PAD",
            "CKD effect with PAD minus CKD effect without PAD"), "exploratory_support")

# GEO explicitly documents already normalized log2 array expression.
d <- load_data("GSE236927", probe = TRUE)
m <- d$m
m$age <- as.numeric(m$`age (y)`)
m$bmi <- as.numeric(m$`bmi (kg(m2)`)
m$fi10 <- as.numeric(m$`fi score (range 0-1)`) / 0.1
m$monocytes <- as.numeric(m$`monocytes (absolute values x109/l)`)
m$lymphocytes <- as.numeric(m$`lymphocytes (absolute values x109/l)`)
sel <- complete.cases(m[, c("age", "bmi", "fi10", "monocytes", "lymphocytes")])
stopifnot(sum(sel) == 88L)
m <- m[sel, ]
for (field in c("age", "bmi", "monocytes", "lymphocytes")) m[[field]] <- m[[field]] - mean(m[[field]])
for (cells in c(FALSE, TRUE)) {
  design <- if (cells) model.matrix(~ fi10 + age + bmi + monocytes + lymphocytes, m) else model.matrix(~ fi10 + age + bmi, m)
  rownames(design) <- m$matrix_column
  name <- if (cells) "GSE236927_fi_cell_adjusted" else "GSE236927_fi"
  e <- d$x[, sel, drop = FALSE]
  keep <- apply(e, 1, var) > 0
  model <- list(e = e[keep, ], a = d$a[keep, ], design = design, m = m,
                block = NULL, correlation = NULL, assay = "processed_log2_array_limma_trend",
                trend = TRUE, input_features = nrow(d$x))
  record_qc(model$e, m, name)
  fit_model(model, single_contrast(design, "fi10", name), name,
            if (cells) "Per +0.1 FI; age/BMI/monocytes/lymphocytes adjusted" else "Per +0.1 FI; age/BMI adjusted",
            "conditional_blood_support")
}

# Renal progression is a different outcome from disability recovery.
d <- load_data("GSE339456")
m <- d$m
m$group <- factor(sub("_.*$", "", m$matrix_column), levels = c("ES", "EP", "LS", "LP"))
stopifnot(!anyNA(m$group), identical(as.integer(table(m$group)), c(10L, 10L, 7L, 8L)))
design <- model.matrix(~ 0 + group, m)
colnames(design) <- c("ES", "EP", "LS", "LP")
rownames(design) <- m$matrix_column
model <- count_input(d, rep(TRUE, nrow(m)), design, "GSE339456_renal_progression")
model$m <- m
contrasts <- makeContrasts(GSE339456_progression_average = ((EP - ES) + (LP - LS)) / 2,
                         GSE339456_early_progression = EP - ES,
                         GSE339456_late_progression = LP - LS,
                         GSE339456_stage_interaction = (LP - LS) - (EP - ES), levels = design)
fit_model(model, contrasts, "GSE339456_renal_progression",
          c("Equal-weight average progression effect across early/late stages",
            "Early progressive minus early stable", "Late progressive minus late stable",
            "Late progression effect minus early progression effect"), "renal_organ_support")

secondary_names <- paste0("GSE144304_worse_", c("sppb", "meangrip", "walktime"))
secondary_tables <- lapply(secondary_names, function(name) read.csv(file.path(out, "tables", paste0(name, "_de.csv"))))
all_secondary <- do.call(rbind, secondary_tables)
all_secondary$FDR_across_three_phenotypes <- p.adjust(all_secondary$P.Value, "BH")
write.csv(all_secondary, file.path(out, "tables", "GSE144304_secondary_joint_fdr.csv"), row.names = FALSE)
all_pathways <- do.call(rbind, lapply(secondary_names, function(name) read.csv(file.path(out, "tables", paste0(name, "_camera.csv")))))
all_pathways$FDR_across_three_phenotypes <- p.adjust(all_pathways$PValue, "BH")
write.csv(all_pathways, file.path(out, "tables", "GSE144304_secondary_pathway_joint_fdr.csv"), row.names = FALSE)
write_json(input_audits, file.path(out, "input_annotation_audit.json"), pretty = TRUE, auto_unbox = TRUE, na = "null")
capture.output(sessionInfo(), file = file.path(out, "analysis_sessionInfo.txt"))
cat("ALL_MODELS_COMPLETE\n")
sink()
