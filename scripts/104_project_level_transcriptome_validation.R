Sys.setlocale("LC_CTYPE", "English_United States.utf8")
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
root <- normalizePath(args[1], winslash = "/", mustWork = TRUE)

suppressPackageStartupMessages({library(dplyr); library(metafor); library(jsonlite)})
options(stringsAsFactors = FALSE)

base <- file.path(root, "muscle_oxphos_reaudit_20261007")
out <- file.path(base, "project_level_validation_20261007")
dir.create(file.path(out, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out, "audit"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out, "reports"), recursive = TRUE, showWarnings = FALSE)

score_path <- file.path(base, "tables", "oxphos_sample_score_results.csv")
gene_path <- file.path(base, "tables", "oxphos_gene_effects_primary6.csv")
score <- read.csv(score_path, stringsAsFactors = FALSE)
genes <- read.csv(gene_path, stringsAsFactors = FALSE)
primary_comparisons <- c("GSE144304_frail_vs_fit", "GSE111006_sarcopenia", "GSE111010_sarcopenia",
                         "GSE111016_sarcopenia", "GSE117525_frailty_baseline", "GSE174106_alm_loss")
family_map <- data.frame(
  accession = c("GSE144304", "GSE111006", "GSE111010", "GSE111016", "GSE117525", "GSE174106"),
  study_family = c("GSE144304", "GSE111017_multicountry", "GSE111017_multicountry",
                   "GSE111017_multicountry", "GSE117525", "GSE174106"),
  stringsAsFactors = FALSE
)

score <- score %>%
  filter(comparison %in% primary_comparisons, variant == "candidate_excluded",
         pathway == "Oxidative Phosphorylation", role == "aging_primary") %>%
  select(-any_of("study_family")) %>%
  left_join(family_map, by = "accession") %>%
  mutate(beta = as.numeric(beta), se = as.numeric(se),
         ci_low = as.numeric(ci_low), ci_high = as.numeric(ci_high),
         p = as.numeric(p), weight = 1 / se^2)
stopifnot(nrow(score) == 6L, all(is.finite(score$beta)), all(is.finite(score$se) & score$se > 0),
          n_distinct(score$study_family) == 4L)

aggregate_family <- function(d) {
  if (nrow(d) == 1L) {
    b <- d$beta; se <- d$se; i2 <- 0; tau2 <- 0; method <- "single_accession"
    low <- d$ci_low; high <- d$ci_high
  } else {
    fit <- rma.uni(yi = d$beta, sei = d$se, method = "REML", test = "knha")
    b <- as.numeric(fit$b); se <- fit$se; i2 <- fit$I2; tau2 <- fit$tau2; method <- "REML_within_project"
    low <- fit$ci.lb; high <- fit$ci.ub
  }
  data.frame(study_family = d$study_family[1], k_accessions = nrow(d), beta = b, se = se,
             ci_low = low, ci_high = high,
             I2_percent = i2, tau2 = tau2, aggregation = method, stringsAsFactors = FALSE)
}

family_effects <- bind_rows(lapply(split(score, score$study_family), aggregate_family))
write.csv(family_effects, file.path(out, "tables", "project_family_effects.csv"), row.names = FALSE)

overall_fit <- rma.uni(yi = family_effects$beta, sei = family_effects$se, method = "REML", test = "knha")
overall <- data.frame(
  analysis = "all_four_project_families", k = overall_fit$k, beta = as.numeric(overall_fit$b), se = overall_fit$se,
  ci_low = overall_fit$ci.lb, ci_high = overall_fit$ci.ub, p = overall_fit$pval,
  I2_percent = overall_fit$I2, tau2 = overall_fit$tau2, stringsAsFactors = FALSE)
write.csv(overall, file.path(out, "tables", "project_family_overall_meta.csv"), row.names = FALSE)

loo <- bind_rows(lapply(family_effects$study_family, function(omit) {
  d <- family_effects %>% filter(study_family != omit)
  fit <- rma.uni(yi = d$beta, sei = d$se, method = "REML", test = "knha")
  data.frame(omitted_project_family = omit, k = fit$k, beta = as.numeric(fit$b), se = fit$se,
             ci_low = fit$ci.lb, ci_high = fit$ci.ub, p = fit$pval, I2_percent = fit$I2,
             tau2 = fit$tau2, stringsAsFactors = FALSE)
}))
write.csv(loo, file.path(out, "tables", "project_family_leave_one_out.csv"), row.names = FALSE)

# Candidate direction stability is summarized at the project-family level. This is a
# validation audit, not a new feature-selection procedure: candidates were fixed before
# this analysis from the six-dataset stability audit.
candidate_genes <- c("COX5A", "GOT2", "MDH2", "SDHA")
genes <- genes %>% filter(comparison %in% primary_comparisons, symbol %in% candidate_genes) %>%
  mutate(accession = sub("_.*", "", comparison)) %>%
  left_join(family_map, by = c("accession" = "accession")) %>%
  mutate(logFC = as.numeric(logFC), p = as.numeric(p), down = logFC < 0)
candidate_family <- genes %>% group_by(symbol, study_family) %>% summarise(
  n_accessions = n(), mean_logFC = mean(logFC, na.rm = TRUE),
  down_fraction = mean(down, na.rm = TRUE), all_down = all(down),
  .groups = "drop")
candidate_summary <- candidate_family %>% group_by(symbol) %>% summarise(
  project_families = n(), families_down = sum(mean_logFC < 0),
  family_down_fraction = mean(mean_logFC < 0),
  families_all_down = sum(all_down), median_family_logFC = median(mean_logFC),
  stable_direction_all_families = all(mean_logFC < 0), .groups = "drop")
candidate_summary$candidate_tier <- ifelse(candidate_summary$stable_direction_all_families,
  "core_directionally_stable", "secondary_directionally_inconsistent")
write.csv(candidate_family, file.path(out, "tables", "candidate_family_direction.csv"), row.names = FALSE)
write.csv(candidate_summary, file.path(out, "tables", "candidate_family_direction_summary.csv"), row.names = FALSE)

manifest <- list(
  generated_utc = format(Sys.time(), tz = "UTC"),
  input_score = score_path, input_gene_effects = gene_path,
  primary_comparisons = primary_comparisons,
  project_families = family_map,
  candidate_genes_fixed_before_validation = candidate_genes,
  score_variant = "candidate_excluded",
  family_aggregation = "REML within family when multiple accessions; single accession retained with reported SE",
  outer_meta = "REML random effects with Knapp-Hartung across four project families",
  interpretation = "Project-family leave-one-out sensitivity, not prospective external validation; no feature selection was redone",
  leakage_control = "GSE111006/GSE111010/GSE111016 treated as one GSE111017 multicountry family"
)
write_json(manifest, file.path(out, "audit", "project_level_validation_specification.json"), pretty = TRUE, auto_unbox = TRUE)

fmt <- function(x, d = 3) ifelse(is.na(x), "NA", formatC(x, format = "f", digits = d))
report <- c(
  "# Project-level transcriptome validation audit",
  "",
  "## Design",
  "",
  "The six primary muscle comparisons were regrouped into four project families. GSE111006, GSE111010 and GSE111016 were treated as one GSE111017 multicountry family because they arise from the same broader project. The OXPHOS candidate-excluded sample-level score and four fixed candidate genes were evaluated without reselecting features. Single-accession confidence intervals retain the original model intervals; the multi-accession family uses its fitted model interval.",
  "",
  "This is a project-family leave-one-out validation sensitivity analysis. It is not prospective external validation because the candidates and score definition were specified in the preceding audit.",
  "",
  "## Family-level effects",
  "",
  "| Project family | Accessions | Beta | 95% CI | I2 (%) | Aggregation |",
  "|---|---:|---:|---:|---:|---|",
  vapply(seq_len(nrow(family_effects)), function(i) paste0("| ", family_effects$study_family[i], " | ", family_effects$k_accessions[i], " | ", fmt(family_effects$beta[i]), " | [", fmt(family_effects$ci_low[i]), ", ", fmt(family_effects$ci_high[i]), "] | ", fmt(family_effects$I2_percent[i], 1), " | ", family_effects$aggregation[i], " |"), character(1)),
  "",
  "## Outer meta-analysis",
  "",
  paste0("Across four project families, beta=", fmt(overall$beta), ", 95% CI ", fmt(overall$ci_low), " to ", fmt(overall$ci_high), ", P=", fmt(overall$p, 4), ", I2=", fmt(overall$I2_percent, 1), "%."),
  "",
  "## Leave-one-project-out analysis",
  "",
  "| Omitted family | K | Beta | 95% CI | P | I2 (%) |",
  "|---|---:|---:|---:|---:|---:|",
  vapply(seq_len(nrow(loo)), function(i) paste0("| ", loo$omitted_project_family[i], " | ", loo$k[i], " | ", fmt(loo$beta[i]), " | [", fmt(loo$ci_low[i]), ", ", fmt(loo$ci_high[i]), "] | ", fmt(loo$p[i], 4), " | ", fmt(loo$I2_percent[i], 1), " |"), character(1)),
  "",
  "The project-family analysis addresses the main independence concern. Because only four families are available, the held-out estimates are precision-limited and should be presented as robustness evidence rather than as a new discovery-validation claim.",
  "",
  "## Fixed candidate direction by project family",
  "",
  "| Gene | Project families | Families with downward mean effect | Downward fraction | All families downward | Median family logFC |",
  "|---|---:|---:|---:|---:|---:|",
  vapply(seq_len(nrow(candidate_summary)), function(i) paste0("| ", candidate_summary$symbol[i], " | ", candidate_summary$project_families[i], " | ", candidate_summary$families_down[i], " | ", fmt(candidate_summary$family_down_fraction[i], 2), " | ", candidate_summary$stable_direction_all_families[i], " | ", fmt(candidate_summary$median_family_logFC[i]), " |"), character(1)),
  "",
  "## Interpretation boundary",
  "",
  "The OXPHOS score remains directionally consistent after collapsing related accessions into project families, but the held-out intervals are wide. The analysis supports reproducibility of a tissue-level association; it does not establish patient-level mediation, causal regulation, or clinical prediction performance."
)
writeLines(report, file.path(out, "reports", "project_level_transcriptome_validation_zh.md"), useBytes = TRUE)
capture.output(sessionInfo(), file = file.path(out, "audit", "R_sessionInfo.txt"))
cat("PROJECT_LEVEL_VALIDATION_COMPLETE\n")
print(overall)
print(loo)
print(candidate_summary)
