Sys.setlocale("LC_CTYPE", "English_United States.utf8")
args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) normalizePath(args[1], winslash = "/") else getwd()

suppressPackageStartupMessages({library(dplyr); library(jsonlite)})
out <- file.path(root, "muscle_oxphos_reaudit_20261007")
dir.create(file.path(out, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out, "audit"), recursive = TRUE, showWarnings = FALSE)

membership <- read.csv(file.path(out, "audit", "oxphos_gene_membership.csv"), stringsAsFactors = FALSE)
primary <- c("GSE144304_frail_vs_fit", "GSE111006_sarcopenia", "GSE111010_sarcopenia",
             "GSE111016_sarcopenia", "GSE117525_frailty_baseline", "GSE174106_alm_loss")
membership <- membership %>% filter(comparison %in% primary)
full <- membership %>% filter(variant == "full")
excluded <- membership %>% filter(variant == "candidate_excluded")
universe <- sort(unique(toupper(full$symbol)))
candidate_genes <- c("COX5A", "GOT2", "MDH2", "SDHA", "CKMT2", "TRIM72")

spec <- data.frame(
  comparison = primary,
  model_path = c(
    file.path("transcriptome_analysis", "models", "GSE144304_older.rds"),
    file.path("transcriptome_expansion_20261006", "models", "GSE111006_sarcopenia.rds"),
    file.path("transcriptome_expansion_20261006", "models", "GSE111010_sarcopenia.rds"),
    file.path("transcriptome_expansion_20261006", "models", "GSE111016_sarcopenia.rds"),
    file.path("transcriptome_expansion_20261006", "models", "GSE117525_frailty_baseline.rds"),
    file.path("transcriptome_expansion_20261006", "models", "GSE174106_alm_loss.rds")
  ), stringsAsFactors = FALSE
)

extract_effects <- function(row) {
  key <- sub("\\.rds$", "", basename(row$model_path))
  x <- readRDS(file.path(root,"input/bulk_expression_inputs.rds"))[[key]]
  e <- x$e; a <- x$a; coef_name <- "contrast"
  fit <- limma::lmFit(e, x$design, weights=x$weights)
  fit <- limma::contrasts.fit(fit, matrix(x$contrast, ncol=1,
    dimnames=list(names(x$contrast),coef_name)))
  fit <- limma::eBayes(fit, trend=x$trend, robust=TRUE)
  a$symbol <- toupper(as.character(a$symbol))
  ok <- which(!is.na(a$symbol) & nzchar(a$symbol) & a$symbol %in% universe)
  ave <- if (!is.null(fit$Amean)) fit$Amean else rowMeans(e)
  ord <- ok[order(-ave[ok], as.character(a$feature_id[ok]))]
  ord <- ord[!duplicated(a$symbol[ord])]
  effect <- fit$coefficients[ord, coef_name]
  pval <- fit$p.value[ord, coef_name]
  data.frame(comparison = row$comparison, symbol = a$symbol[ord], logFC = as.numeric(effect),
             p = as.numeric(pval), stringsAsFactors = FALSE)
}

effects <- bind_rows(lapply(seq_len(nrow(spec)), function(i) extract_effects(spec[i, ])))
effects <- effects %>% group_by(comparison) %>% mutate(fdr = p.adjust(p, method = "BH")) %>% ungroup() %>%
  mutate(direction_down = logFC < 0)

membership_summary <- lapply(universe, function(g) {
  f <- full %>% filter(toupper(symbol) == g)
  e <- excluded %>% filter(toupper(symbol) == g)
  z <- effects %>% filter(symbol == g)
  data.frame(
    symbol = g, candidate_gene = g %in% candidate_genes,
    member_n_full = n_distinct(f$comparison), member_n_candidate_excluded = n_distinct(e$comparison),
    covered_n_effect = n_distinct(z$comparison), down_n = sum(z$direction_down, na.rm = TRUE),
    down_fraction = if (nrow(z)) mean(z$direction_down) else NA_real_,
    median_logFC = if (nrow(z)) median(z$logFC, na.rm = TRUE) else NA_real_,
    min_logFC = if (nrow(z)) min(z$logFC, na.rm = TRUE) else NA_real_,
    max_logFC = if (nrow(z)) max(z$logFC, na.rm = TRUE) else NA_real_,
    nominal_p05_n = sum(z$p < .05, na.rm = TRUE), fdr05_n = sum(z$fdr < .05, na.rm = TRUE),
    stable_member = n_distinct(e$comparison) >= 5,
    stable_direction = nrow(z) >= 5 && mean(z$direction_down) >= 2/3,
    stable_core_member = n_distinct(e$comparison) >= 5 && nrow(z) >= 5 && mean(z$direction_down) >= 2/3,
    stringsAsFactors = FALSE
  )
}) %>% bind_rows() %>% arrange(desc(stable_core_member), desc(member_n_candidate_excluded), desc(down_fraction))

write.csv(effects, file.path(out, "tables", "oxphos_gene_effects_primary6.csv"), row.names = FALSE)
write.csv(membership_summary, file.path(out, "tables", "oxphos_gene_member_stability.csv"), row.names = FALSE)
write.csv(membership_summary %>% filter(candidate_gene), file.path(out, "tables", "oxphos_candidate_gene_stability.csv"), row.names = FALSE)

stable <- membership_summary %>% filter(stable_core_member)
summary <- list(
  primary_datasets = primary, gene_universe_n = length(universe), stable_core_n = nrow(stable),
  stable_core_genes = stable$symbol, candidate_genes = candidate_genes,
  candidate_excluded_before_effect_analysis = TRUE,
  direction_definition = "logFC < 0 in the adverse phenotype contrast",
  stable_core_definition = "present after candidate exclusion in >=5/6 primary datasets and downward in >=2/3 covered datasets",
  interpretation = "gene-level directions are descriptive across heterogeneous studies; the sample-level module score remains the primary inferential result"
)
write_json(summary, file.path(out, "audit", "oxphos_gene_stability_summary.json"), pretty = TRUE, auto_unbox = TRUE)
writeLines(capture.output(print(head(membership_summary, 25))), file.path(out, "audit", "oxphos_gene_stability_top25.txt"), useBytes = TRUE)
cat("OXPHOS_GENE_STABILITY_COMPLETE\n")
print(head(membership_summary, 20))
