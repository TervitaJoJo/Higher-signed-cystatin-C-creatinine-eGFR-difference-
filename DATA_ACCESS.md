# Data access and submission boundary

The local package contains processed observations for complete model reproduction. These are real processed data, not synthetic substitutes. Clinical data were derived from HRS 2016 VBS, RAND HRS 1992-2022 v1, the 2022 Tracker, CHARLS 2011 blood data and Harmonized CHARLS version D. Public molecular data retain GEO accessions for scientific provenance.

HRS policy requires users to obtain HRS data through the authorised access process. Original survey files and derived participant-level extracts must not be placed in a public journal supplement or transferred to reviewers merely because identifiers were relabelled. The HRS replication policy explicitly states that investigators should not provide HRS data to journals for redistribution. Relevant official policies:

- https://hrs.isr.umich.edu/data-products/data-access/conditions-of-use
- https://hrs.isr.umich.edu/sites/default/files/HRSReplicationRules_v5b.pdf
- https://charls.charlsdata.com/
- https://www.g2aging.org/

The `local_data` directory is therefore excluded by `make_submission.py`. The submission archive contains executable analysis code, clinical aggregate outputs and public processed molecular observations. It supports figure recreation and public molecular model refitting. Complete clinical model refitting requires the recipient's authorised survey access and the processed-input contract described below. This distinction is also reported in `run_summary.json`.

## Processed clinical input contract

`clinical_panels.rds`: named R list `CHARLS2011`, `HRS2016VBS`; each contains `baseline`, `panel` and the original wave-range metadata. The panel retains the original analysis variables and observed wave states before next-wave risk-set filtering. `pid` is a within-cohort project index. PSU and stratum labels are recoded factors with their original level ordering.

`renal_intervals.rds`: named list with the same cohorts; exact final primary intervals augmented with the already linked BUN, haemoglobin and HRS serum albumin. All numerical marker values and missingness are retained. Clinical model scripts reconstruct covariates and fit survey-weighted models from these observations.

`selection_intervals.rds`: the original known-outcome, creatinine-eligible interval data used for the complete-case selection sensitivity analysis. This wider denominator is distinct from the final signed-difference risk set.

Raw HHID, PN, ID, wave-1 ID, linking keys, labels and source-file metadata are excluded. No linking key or reverse lookup is included. The files remain pseudonymised survey microdata and are for authorised local use. Do not submit the entire local `code` directory; submit the archive produced by `make_submission.py`.

## Public processed molecular input contract

`bulk_expression_inputs.rds`: 13 named analysis inputs. Each supplies gene-by-sample expression, any voom weights, feature annotations, fixed design and contrast, and public series provenance. These are observations/designs, not archived fitted coefficients. The six primary contrasts are grouped into four project families by the analysis code. The inputs also retain the existing secondary and sensitivity contrasts used by the OXPHOS audit.

`muscle_pseudobulk_counts.npz`: public atlas counts aggregated before the present donor-level analysis. Rows match `muscle_pseudobulk_metadata.csv` exactly. Metadata retain only anonymous donor codes, cell-group annotation, age, sex and cell count. The atlas source is CELLxGENE dataset `b3c9904b-2ffd-4038-8620-ee1c5b072e9c`.

`figure_source_data.json`: CSV-formatted aggregate tables containing the risk-set flow, exposure quantiles, exact histogram bin counts and archived spline curves. It contains no person identifiers or row-level clinical exposure list.

`reference_tables.json`: complete archived numeric result tables, used as fixed comparison targets and as aggregate drawing inputs. It is never evidence of model refitting by itself. The runner removes targets that must be recomputed before running models, so missing computations cannot pass by reusing a seeded target.
