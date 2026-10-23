Sys.setlocale("LC_CTYPE", "English_United States.utf8")
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
root <- normalizePath(args[1], winslash = "/", mustWork = TRUE)

suppressPackageStartupMessages({library(limma); library(jsonlite); library(metafor)})
options(stringsAsFactors = FALSE)
set.seed(20261007)
out <- file.path(root, "muscle_oxphos_reaudit_20261007")
for (folder in c("tables", "audit", "data_internal", "figures")) {
  dir.create(file.path(out, folder), recursive = TRUE, showWarnings = FALSE)
}
sink(file.path(out, "audit", "analysis_run.log"), split = TRUE)
excluded <- c("COX5A", "GOT2", "MDH2", "SDHA", "CKMT2", "TRIM72")
ox <- "Oxidative Phosphorylation"
lines <- strsplit(readLines(file.path(root, "transcriptome_analysis", "resources",
                                     "MSigDB_Hallmark_2020.gmt")), "\t")
sets <- lapply(lines, function(x) unique(x[-c(1, 2)]))
names(sets) <- vapply(lines, `[`, character(1), 1L)
stopifnot(length(sets) == 50L, ox %in% names(sets))
spec <- data.frame(
  comparison = c("GSE144304_frail_vs_fit", "GSE144304_worse_sppb", "GSE144304_worse_meangrip",
                 "GSE144304_worse_walktime", "GSE111006_sarcopenia", "GSE111010_sarcopenia",
                 "GSE111016_sarcopenia", "GSE117525_frailty_baseline", "GSE174106_alm_loss",
                 "GSE117525_frailty_title_prefix", "GSE157712_ckd_vs_control",
                 "GSE287118_ckd_no_pad", "GSE157712_ckd_race_adjusted"),
  source = c(rep("transcriptome_analysis", 4), rep("transcriptome_expansion_20261006", 6),
             rep("transcriptome_analysis", 3)),
  model = c("GSE144304_older", "GSE144304_worse_sppb", "GSE144304_worse_meangrip",
            "GSE144304_worse_walktime", "GSE111006_sarcopenia", "GSE111010_sarcopenia",
            "GSE111016_sarcopenia", "GSE117525_frailty_baseline", "GSE174106_alm_loss",
            "GSE117525_frailty_title_prefix", "GSE157712_independent",
            "GSE287118_no_pad", "GSE157712_race_sensitivity"),
  role = c("aging_primary", rep("functional_secondary", 3), rep("aging_primary", 5),
           "batch_proxy_sensitivity", "renal_primary", "renal_processed_support",
           "renal_covariate_sensitivity")
)
spec$accession <- sub("_.*", "", spec$comparison)
spec$study_family <- ifelse(spec$accession %in% c("GSE111006", "GSE111010", "GSE111016"),
                          "GSE111017_multicountry", spec$accession)
spec$original_nine <- seq_len(nrow(spec)) <= 9L
write.csv(spec, file.path(out, "audit", "analysis_manifest.csv"), row.names = FALSE)
write_json(list(status = "exploratory_reaudit_after_prior_results", date = "2026-10-07",
  pathway_library = "MSigDB_Hallmark_2020", primary_variant = "candidate_excluded",
  methods = c("CAMERA_estimated_correlation_competitive", "FRY_directional_self_contained",
              "CAMERA_fixed_0.01_sensitivity", "sample_score_HC3"),
  score = "Mean within-cohort gene z-scores, standardized to one sample SD",
  FDR = "Each method/variant: within-comparison Hallmark50, all primary6 x Hallmark50, and all original9 x Hallmark50",
  meta = "Exploratory heterogeneous adverse-phenotype group contrasts; REML with Knapp-Hartung; BH over Hallmark50",
  excluded_candidates = excluded, no_patient_level_linkage = TRUE),
  file.path(out, "audit", "analysis_specification.json"), pretty = TRUE, auto_unbox = TRUE)

load_model <- function(i) {
  input <- readRDS(file.path(root, "input/bulk_expression_inputs.rds"))[[spec$model[i]]]
  d <- input
  e <- input$e; weights <- input$weights; contrast <- input$contrast; trend <- input$trend
  stopifnot(identical(colnames(e), as.character(d$m$matrix_column)),
    identical(colnames(e), rownames(d$design)), all(is.finite(e)))
  a <- d$a
  a$symbol <- toupper(as.character(a$symbol))
  idx <- which(!is.na(a$symbol) & nzchar(a$symbol))
  idx <- idx[order(-rowMeans(e)[idx], a$feature_id[idx])]
  idx <- idx[!duplicated(a$symbol[idx])]
  e <- e[idx, , drop = FALSE]
  rownames(e) <- a$symbol[idx]
  if (!is.null(weights)) weights <- weights[idx, , drop = FALSE]
  list(e = e, weights = weights, design = d$design, contrast = contrast,
       m = d$m, trend = trend, feature_id = a$feature_id[idx])
}

# Scores retain samples, rather than genes, as the independent units of inference.
score_fit <- function(y, design, contrast) {
  fit <- lm.fit(design, y)
  stopifnot(fit$rank == ncol(design))
  bread <- solve(crossprod(design))
  leverage <- rowSums((design %*% bread) * design)
  adjusted_residuals <- fit$residuals / (1 - leverage)
  meat <- crossprod(design, design * adjusted_residuals^2)
  vcov <- bread %*% meat %*% bread
  beta <- sum(contrast * fit$coefficients)
  se <- sqrt(drop(t(contrast) %*% vcov %*% contrast))
  df <- nrow(design) - ncol(design)
  p <- 2 * pt(-abs(beta / se), df)
  bound <- qt(0.975, df) * se
  list(beta = beta, se = se, ci_low = beta - bound, ci_high = beta + bound,
       p = p, df = df, max_leverage = max(leverage),
       max_abs_studentized = max(abs(fit$residuals / sqrt(sum(fit$residuals^2) / df * (1 - leverage)))))
}
score_vector <- function(e, genes, min_genes = 15L) {
  idx <- which(rownames(e) %in% genes & apply(e, 1L, sd) > 1e-8)
  if (length(idx) < min_genes) return(NULL)
  x <- e[idx, , drop = FALSE]
  z <- (x - rowMeans(x)) / apply(x, 1L, sd)
  raw <- colMeans(z)
  list(score = as.numeric(scale(raw)), genes = rownames(x), raw = raw)
}
models <- lapply(seq_len(nrow(spec)), load_model)
primary_idx <- which(spec$role == "aging_primary")
common_ox <- Reduce(intersect, lapply(models[primary_idx], function(d) {
  intersect(rownames(d$e), setdiff(sets[[ox]], excluded))
}))
write.csv(data.frame(symbol = common_ox), file.path(out, "audit", "common_oxphos_genes.csv"), row.names = FALSE)
overlap <- list()
for (i in seq_len(nrow(spec) - 1L)) for (j in seq.int(i + 1L, nrow(spec))) {
  n_shared <- length(intersect(models[[i]]$m$sample_accession, models[[j]]$m$sample_accession))
  overlap[[length(overlap) + 1L]] <- data.frame(comparison_a = spec$comparison[i],
    comparison_b = spec$comparison[j], shared_GSM = n_shared)
}
write.csv(do.call(rbind, overlap), file.path(out, "audit", "GSM_sample_overlap.csv"), row.names = FALSE)
path_rows <- score_rows <- sample_rows <- qc_rows <- membership <- proxy_rows <- loo_rows <- list()
for (i in seq_len(nrow(spec))) {
  d <- models[[i]]
  qc_rows[[i]] <- data.frame(comparison = spec$comparison[i], n = ncol(d$e),
    genes = nrow(d$e), rank = qr(d$design)$rank, residual_df = nrow(d$design) - ncol(d$design),
    condition_number = kappa(d$design), group_counts = paste(names(table(d$m$group)), table(d$m$group), collapse = ";"),
    covariates = paste(colnames(d$design), collapse = ";"), original_model_weights = !is.null(d$weights),
    sample_alignment_passed = TRUE)
  for (variant in c("full", "candidate_excluded")) {
    keep <- if (variant == "full") rep(TRUE, nrow(d$e)) else !rownames(d$e) %in% excluded
    e <- d$e[keep, , drop = FALSE]
    weights <- if (is.null(d$weights)) NULL else d$weights[keep, , drop = FALSE]
    used_sets <- if (variant == "full") sets else lapply(sets, setdiff, y = excluded)
    indices <- lapply(used_sets, function(s) which(rownames(e) %in% s))
    indices <- indices[lengths(indices) >= 15L]
    tests <- list(
      CAMERA_estimated = camera(e, indices, d$design, d$contrast, weights = weights,
                               inter.gene.cor = NA, trend.var = d$trend, sort = FALSE),
      FRY_directional = fry(e, indices, d$design, d$contrast, weights = weights,
                            trend.var = d$trend, robust = TRUE, sort = "none"),
      CAMERA_fixed_0.01 = camera(e, indices, d$design, d$contrast, weights = weights,
                                inter.gene.cor = 0.01, trend.var = d$trend, sort = FALSE)
    )
    for (method in names(tests)) {
      tab <- tests[[method]]
      extra <- data.frame(comparison = spec$comparison[i], accession = spec$accession[i],
                          role = spec$role[i], original_nine = spec$original_nine[i], variant = variant,
                          method = method, pathway = rownames(tab), genes = tab$NGenes,
                          direction = tab$Direction, p = tab$PValue,
                          estimated_gene_correlation = if ("Correlation" %in% names(tab)) tab$Correlation else NA_real_,
                          p_mixed = if ("PValue.Mixed" %in% names(tab)) tab$PValue.Mixed else NA_real_)
      extra$fdr_hallmark_within_comparison <- p.adjust(extra$p, "BH")
      path_rows[[length(path_rows) + 1L]] <- extra
    }
    for (pathway in names(indices)) {
      sc <- score_vector(e, used_sets[[pathway]])
      fit <- score_fit(sc$score, d$design, d$contrast)
      score_rows[[length(score_rows) + 1L]] <- cbind(data.frame(comparison = spec$comparison[i],
        accession = spec$accession[i], study_family = spec$study_family[i], role = spec$role[i],
        original_nine = spec$original_nine[i], variant = variant, pathway = pathway,
        n = ncol(e), genes = length(sc$genes)), as.data.frame(fit))
      if (pathway == ox) {
        membership[[length(membership) + 1L]] <- data.frame(comparison = spec$comparison[i], variant = variant,
                                                          symbol = sc$genes)
        sample_rows[[length(sample_rows) + 1L]] <- data.frame(comparison = spec$comparison[i], variant = variant,
          sample_accession = d$m$sample_accession, matrix_column = d$m$matrix_column,
          group = as.character(d$m$group), score_sd = sc$score, mean_gene_z = sc$raw)
      }
    }
  }
  e <- d$e[!rownames(d$e) %in% excluded, , drop = FALSE]
  sc <- score_vector(e, common_ox)
  fit <- score_fit(sc$score, d$design, d$contrast)
  score_rows[[length(score_rows) + 1L]] <- cbind(data.frame(comparison = spec$comparison[i],
    accession = spec$accession[i], study_family = spec$study_family[i], role = spec$role[i],
    original_nine = spec$original_nine[i], variant = "common_gene_sensitivity", pathway = ox,
    n = ncol(e), genes = length(sc$genes)), as.data.frame(fit))
  base <- score_vector(e, setdiff(sets[[ox]], excluded))$score
  proxy_sets <- list(slow_fiber = c("MYH7", "TNNT1", "TNNI1"),
                     fast_fiber = c("MYH1", "MYH2", "TNNT3", "TNNI2"),
                     immune = c("PTPRC", "LST1", "TYROBP"))
  proxies <- lapply(proxy_sets, function(g) score_vector(e, g, min_genes = 2L))
  available <- !vapply(proxies, is.null, logical(1))
  if (all(available)) {
    design_proxy <- cbind(d$design, do.call(cbind, lapply(proxies, `[[`, "score")))
    colnames(design_proxy)[seq.int(ncol(d$design) + 1L, ncol(design_proxy))] <- names(proxies)
    if (qr(design_proxy)$rank == ncol(design_proxy)) {
      fit <- score_fit(base, design_proxy, c(d$contrast, rep(0, length(proxies))))
      proxy_rows[[length(proxy_rows) + 1L]] <- cbind(data.frame(comparison = spec$comparison[i],
        role = spec$role[i], status = "ok", proxy_adjustment = "slow_fast_immune_expression_proxies"), as.data.frame(fit))
    }
  }
  if (spec$role[i] == "aging_primary" || spec$role[i] == "functional_secondary") {
    for (j in seq_len(ncol(e))) {
      x <- d$design[-j, , drop = FALSE]
      if (qr(x)$rank != ncol(x)) next
      fit <- score_fit(base[-j], x, d$contrast)
      loo_rows[[length(loo_rows) + 1L]] <- cbind(data.frame(comparison = spec$comparison[i],
        omitted_GSM = d$m$sample_accession[j]), as.data.frame(fit))
    }
  }
  cat("COMPLETE", spec$comparison[i], "n", ncol(e), "OXPHOS_genes", length(sc$genes), "\n")
}

paths <- do.call(rbind, path_rows)
paths$fdr_primary6_all_hallmark <- paths$fdr_original9_all_hallmark <- NA_real_
for (method in unique(paths$method)) for (variant in unique(paths$variant)) {
  ids <- which(paths$method == method & paths$variant == variant & paths$role == "aging_primary")
  paths$fdr_primary6_all_hallmark[ids] <- p.adjust(paths$p[ids], "BH")
  ids <- which(paths$method == method & paths$variant == variant & paths$original_nine)
  paths$fdr_original9_all_hallmark[ids] <- p.adjust(paths$p[ids], "BH")
}
scores <- do.call(rbind, score_rows)
scores$fdr_hallmark_within_comparison <- scores$fdr_primary6_all_hallmark <- scores$fdr_original9_all_hallmark <- NA_real_
for (variant in unique(scores$variant)) {
  for (name in unique(scores$comparison)) {
    ids <- which(scores$variant == variant & scores$comparison == name)
    scores$fdr_hallmark_within_comparison[ids] <- p.adjust(scores$p[ids], "BH")
  }
  ids <- which(scores$variant == variant & scores$role == "aging_primary")
  scores$fdr_primary6_all_hallmark[ids] <- p.adjust(scores$p[ids], "BH")
  ids <- which(scores$variant == variant & scores$original_nine)
  scores$fdr_original9_all_hallmark[ids] <- p.adjust(scores$p[ids], "BH")
}
write.csv(paths, file.path(out, "tables", "hallmark_correlation_aware_results.csv"), row.names = FALSE)
write.csv(subset(paths, pathway == ox), file.path(out, "tables", "oxphos_correlation_aware_results.csv"), row.names = FALSE)
write.csv(scores, file.path(out, "tables", "hallmark_sample_score_results.csv"), row.names = FALSE)
write.csv(subset(scores, pathway == ox), file.path(out, "tables", "oxphos_sample_score_results.csv"), row.names = FALSE)
write.csv(do.call(rbind, sample_rows), file.path(out, "data_internal", "oxphos_public_sample_scores.csv"), row.names = FALSE)
write.csv(do.call(rbind, membership), file.path(out, "audit", "oxphos_gene_membership.csv"), row.names = FALSE)
write.csv(do.call(rbind, qc_rows), file.path(out, "audit", "model_alignment_and_design.csv"), row.names = FALSE)
proxies <- do.call(rbind, proxy_rows)
proxies$fdr_all_proxy_tests <- p.adjust(proxies$p, "BH")
write.csv(proxies, file.path(out, "tables", "oxphos_composition_proxy_sensitivity.csv"), row.names = FALSE)
write.csv(do.call(rbind, loo_rows), file.path(out, "tables", "oxphos_leave_one_sample_out.csv"), row.names = FALSE)

meta_rows <- meta_loo <- list()
meta_fit <- function(d, variant, pathway, analysis, omitted = "") {
  if (nrow(d) < 3L) return(NULL)
  fit <- rma.uni(yi = d$beta, sei = d$se, method = "REML", test = "knha")
  data.frame(variant = variant, pathway = pathway, analysis = analysis, omitted = omitted,
             k = fit$k, beta = as.numeric(fit$b), se = fit$se, ci_low = fit$ci.lb, ci_high = fit$ci.ub,
             p = fit$pval, tau2 = fit$tau2, I2 = fit$I2)
}
for (variant_value in unique(scores$variant)) for (pathway_value in unique(scores$pathway)) {
  d <- scores[scores$role == "aging_primary" & scores$variant == variant_value &
              scores$pathway == pathway_value, , drop = FALSE]
  if (nrow(d) < 3L) next
  meta_rows[[length(meta_rows) + 1L]] <- meta_fit(d, variant_value, pathway_value, "one_contrast_per_cohort")
  # Grouping the three country cohorts tests sensitivity to a shared research project.
  grouped <- do.call(rbind, lapply(split(d, d$study_family), function(g) {
    w <- 1 / g$se^2
    data.frame(beta = sum(w * g$beta) / sum(w), se = sqrt(1 / sum(w)))
  }))
  meta_rows[[length(meta_rows) + 1L]] <- meta_fit(grouped, variant_value, pathway_value, "shared_project_grouped_sensitivity")
  if (pathway_value == ox) {
    for (accession in unique(d$accession)) {
      left <- d[d$accession != accession, , drop = FALSE]
      meta_loo[[length(meta_loo) + 1L]] <- meta_fit(left, variant_value,
        pathway_value, "leave_one_cohort_out", accession)
    }
    for (family in unique(d$study_family)) {
      left <- d[d$study_family != family, , drop = FALSE]
      meta_loo[[length(meta_loo) + 1L]] <- meta_fit(left, variant_value,
        pathway_value, "leave_one_project_out", family)
    }
  }
}
meta <- do.call(rbind, meta_rows)
meta$fdr_all_hallmark <- NA_real_
for (variant in unique(meta$variant)) for (analysis in unique(meta$analysis)) {
  ids <- which(meta$variant == variant & meta$analysis == analysis)
  meta$fdr_all_hallmark[ids] <- p.adjust(meta$p[ids], "BH")
}
write.csv(meta, file.path(out, "tables", "sample_score_exploratory_meta.csv"), row.names = FALSE)
write.csv(do.call(rbind, meta_loo), file.path(out, "tables", "oxphos_meta_leave_one_out.csv"), row.names = FALSE)
write_json(list(completed = TRUE, comparisons = nrow(spec), primary_cohorts = length(primary_idx),
  unique_primary_GSM = length(unique(unlist(lapply(models[primary_idx], function(d) d$m$sample_accession)))),
  sum_primary_samples = sum(vapply(models[primary_idx], function(d) ncol(d$e), integer(1))),
  shared_research_project = "GSE111006/GSE111010/GSE111016 are GSE111017 subseries with PMID31862890",
  common_oxphos_genes = length(common_ox), gene_set_tests = nrow(paths), score_tests = nrow(scores),
  source_model_checks_passed = TRUE, clinical_discordance_available_in_muscle = FALSE),
  file.path(out, "audit", "run_summary.json"), pretty = TRUE, auto_unbox = TRUE)
capture.output(sessionInfo(), file = file.path(out, "audit", "R_sessionInfo.txt"))
cat("MUSCLE_OXPHOS_REAUDIT_COMPLETE\n")
sink()
