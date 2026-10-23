# Anonymous analysis and reproduction package

This package accompanies the study of next-wave recovery after established activities of daily living difficulty. HRS is the primary cohort; CHARLS is an exploratory comparison. External muscle analyses provide unlinked biological context.

## One-command entry points

Open a terminal in this directory. With Python dependencies installed, double-click `run_review.cmd`, or run:

```powershell
.\run_review.cmd
```

The default `figures` mode recreates the three main figures and six supplementary figures from locked aggregate sources. It does not refit clinical models. The exposure histogram uses exact archived bin counts, not participant-level records. Result numbers are preserved; supplementary figure styling can differ from the manuscript artwork.

```powershell
.\run_review.cmd --mode full
.\run_review.cmd --mode public
```

`full` reconstructs adjacent-wave risk sets from the processed longitudinal panels, refits clinical models and sensitivity analyses, runs the original 300 survey-PSU bootstrap samples for the internal discrimination checks, recomputes external muscle results, compares reproduced tables with locked outputs, and redraws figures. It requires the locally held processed clinical inputs in `local_data`. Double-click `run_full.cmd` for this complete local workflow.

`public` recomputes the external bulk-muscle and donor-level analyses from bundled processed expression data, checks their results, and recreates all figures using the aggregate clinical source data. This mode works in the submission archive, which excludes `local_data`.

Cross-platform entry: `python reproduce.py --mode figures` (or `full`, `public`). Results are written to `results`. A temporary directory outside the package holds intermediate models and logs and is removed after success. Failed runs retain that directory for diagnosis. `--work-dir <directory> --keep-work` preserves a selected external working directory.

## Dependencies

Python package versions used for verification are in `requirements.txt`:

```powershell
python -m pip install -r requirements.txt
Rscript scripts/setup.R
```

The runner automatically uses the default library created by `setup.R`. Set `REVIEW_R_LIB` for a different library. Set `RSCRIPT` to the executable if Rscript is not on PATH. Dependencies are installed outside this package. R 4.6.0 was used for the verification run. `data/original_R_session.txt` records statistical package versions. Computational results are checked with numerical tolerances because BLAS, operating systems and package versions can differ.

## Contents

- `scripts`: portable analysis implementations and plotting functions.
- `raw_pipeline`: original source-data construction, public-expression acquisition, processing and model-building code, with anonymous configurable paths. See `PIPELINE.md` for its execution boundary and source requirements.
- `data`: processed public molecular inputs, clinical aggregate source data, locked reference outputs and dictionary.
- `local_data`: pseudonymised clinical panels and intervals retained for authorised local reproduction; excluded from the submission archive.
- `results`: regenerated figures, comparison results and run summary.
- `make_submission.py`: builds `Anonymous_peer_review_code.zip` in the parent directory and checks file count and common identity markers. The archive and the local source tree are separate, so there is no duplicated archive nested inside the 100-file package.

The package contains no fitted clinical model objects, original survey IDs, author names, email addresses or personal absolute paths. Relabelled participants and PSUs retain repeated-observation and survey-design relationships. Relabelling does not make restricted survey microdata freely redistributable. Follow `DATA_ACCESS.md` when submitting.

All original manuscript and supplementary documents remain outside this package and are unchanged. No new analysis is introduced by packaging. A comparison failure stops the runner before it declares success.
