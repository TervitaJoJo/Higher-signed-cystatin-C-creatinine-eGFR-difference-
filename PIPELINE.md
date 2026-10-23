# Analysis provenance and execution map

## Executed processed-data pipeline

1. `prepare_clinical.R`: materialise original panel interfaces; construct adjacent-wave intervals and primary complete cases using the original functions.
2. `102_clinical_nonlinear_interaction_competing_risk.R`: survey-weighted quasibinomial models for recovery, persistent difficulty and death; original spline basis, interaction screen, standardised probability contrasts and exploratory synthesis.
3. `materialise_primary.R`: format the six original state-model coefficients for the main figure source table. This is a representation change, not an additional model.
4. `73_complete_case_selection.R`: original complete-case selection weights and sensitivity fits, using the original wider known-outcome denominator.
5. `106_mortality_heterogeneity.R`, `113_hrs_primary_evalue_and_framing.R`: existing secondary state and sensitivity calculations.
6. `119_renal_strata_discordance.R`, `124_renal_extension.R`: original filtration strata, relative categories, joint phenotypes, combined filtration and matched blood-marker adjustments. The processed-data adapter replaces raw ID linkage with the exact previously linked marker observations.
7. `130_primary_baseline_characteristics.R`: unique-person descriptive characteristics.
8. `105_incremental_prediction_gain.R`, `110_prediction_gain_design_validation.R`: first-eligible-interval discrimination checks, including the original PSU-grouped folds and 300 stratified rescaled bootstrap samples. The former provides the original first-interval cache; the latter produces the survey-design-validated check reported in the supplement.
9. `export_clinical_visuals.R`: refit the original recovery spline and check its coefficients, events and nonlinearity tests before exporting pointwise curves.
10. `98_muscle_oxphos_reaudit.R`: original correlation-aware pathway tests, candidate-excluded OXPHOS sample scoring and HC3 sample-level contrasts. `101_muscle_oxphos_gene_stability.R` refits the gene models from the processed expression and fixed designs. `104_project_level_transcriptome_validation.R` recomputes the four-family synthesis and leave-one-family-out checks.
11. `single_cell_oxphos.py`, `single_cell_lr.py`: original donor-level age associations, DoRothEA A/B/C ULM screen and expression co-availability screen from public pseudobulk observations. The DoRothEA network is frozen in `data/dorothea_network.csv` and checked by reproducing the archived regulator estimates.
12. Plot scripts reconstruct Figures 1-3 and S1-S6. Figure 3B uses the current matched-case estimate matrix. The public clinical histogram uses the exact bins and counts from the original interval-level histogram.

## Source-data construction code

The `raw_pipeline` directory preserves the earlier source acquisition, survey harmonisation, ADL coding, death propagation, biomarker equations, expression processing and limma model construction. This layer documents the full provenance before the bundled processed-data starting point. It has not been substituted with synthetic data or a new statistical specification.

Survey construction order: copy the desired raw scripts plus `scripts/renal_equations.R` into an external workspace's `code` directory; create `audit`, `tables`, `figures`, `logs`, `reports` and `data_internal`; set the environment variables below; run `02_build_charls_panel.R`, `07_build_hrs_renal_panel.R`, `08_common_endpoint_renal_models.R`, then `10_marker_discordance_analysis.R`, passing the workspace root as their argument. These historical scripts also construct earlier comparison samples not used in the final manuscript. Their full source-file set is required when executing this historical layer.

- `CHARLS_RAW_DIR`: source root matching the original relative CHARLS release layout.
- `HRS_RAND_FILE`: RAND HRS 1992-2022 v1 Stata file.
- `HRS_TRACKER_FILE`: 2022 Tracker Stata file.
- `HRS_VBS_FILE`: 2016 VBS SPSS file.
- `HRS_DBS_FILE`: 2010 biomarker Stata file required by the historical construction script.
- `CHARLS_BLOOD_FILE`: 2011 blood Stata file, including BUN and haemoglobin.
- `CHARLS_HARMONISED_FILE`: Harmonized CHARLS version D Stata file.

The two raw survey construction scripts also save private linkage caches in `data_internal/private_linkage` for the raw renal-marker script. These caches remain in the authorised raw workspace. The tested one-command pipeline starts from the already linked, pseudonymised observations and does not need those identifiers.

Molecular preparation order: public metadata (`38`), source acquisition and clinical metadata audit (`39`, `40`, `43`), primary expression models (`44`); expanded public metadata, acquisition and preparation (`49`-`52`), then expanded models (`53`). Raw rebuilding also requires the corresponding platform annotations and Bioconductor annotation packages used in those scripts. Network resources can change; the bundled processed matrices and fixed designs are the frozen starting point for exact verification. Atlas pseudobulk construction is in `108_multidonor_muscle_pseudobulk.py` and requires the source H5AD file. Public dataset IDs remain in the code as third-party provenance.

## Preserved interpretation

All clinical models estimate adjacent-wave odds ratios. Death is an observed competing state and is included as non-recovery in recovery models. Arithmetic mean eGFR and combined creatinine-cystatin C eGFR remain separate adjustment specifications. No continuous-time competing-risk or treatment analysis is introduced. CHARLS remains exploratory. Molecular data remain unlinked and do not establish mediation.
