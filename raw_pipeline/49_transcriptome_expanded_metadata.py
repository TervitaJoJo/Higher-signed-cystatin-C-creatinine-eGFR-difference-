"""Audit phenotype-first shortlist before examining new expression results."""
from __future__ import annotations

import csv
import importlib.util
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "transcriptome_expansion_20261006"
RAW = OUT / "source_records"
CANDIDATES = {
    "GSE111006": "Independent UK older-male muscle sarcopenia population; parent not counted again",
    "GSE111010": "Independent Jamaica older-male muscle sarcopenia population; parent not counted again",
    "GSE111016": "Independent Singapore older-male muscle sarcopenia population; parent not counted again",
    "GSE117525": "Larger muscle cohort with age/muscle weakness; verify public function linkage",
    "GSE226151": "Sarcopenia multiomics; verify human tissue, actual independent patients and RNA matrix",
    "GSE174106": "Human muscle aging; verify public strength/function phenotype",
    "GSE140358": "Independent frail/robust peripheral blood; tissue-compatible possible PCR route",
    "GSE287726": "Frailty in older cardiovascular patients; verify covariates and genomic availability",
    "GSE129534": "Race-specific middle-age frailty; verify age/race confounding and cell type",
    "GSE15072": "CKD PBMC; independent renal blood expression support",
    "GSE37171": "Uremic peripheral blood; check batch, independent donors, dialysis and controls",
    "GSE62792": "CKD of unknown etiology whole blood; separate endemic-region/environmental influences",
    "GSE142153": "Diabetic nephropathy/ESRD PBMC; verify glycemic and renal covariate availability",
    "GSE230848": "DKD progression and blood/renal signals; verify actual deposited material",
    "GSE298126": "Older multidomain intervention; verify paired blood expression and function phenotype",
    "GSE40366": "Vitality 90+ blood; verify frailty/disability and donor/timepoint labels",
    "GSE65218": "Vitality 90+ mortality; explicitly check overlap with GSE40366",
}


def load_metadata_helpers():
    spec = importlib.util.spec_from_file_location("expanded_geo_metadata", Path(__file__).with_name("38_transcriptome_candidate_metadata_audit.py"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    module.OUT = OUT
    module.RAW = RAW
    return module


def main() -> None:
    sys.stdout.reconfigure(encoding="utf-8")
    RAW.mkdir(parents=True, exist_ok=True)
    helpers = load_metadata_helpers()
    audit = []
    all_samples = []
    for accession, reason in CANDIDATES.items():
        series_text, series_provenance = helpers.fetch(accession, "self")
        series = helpers.parse_soft(series_text, "SERIES")
        series = series[0] if series else {}
        sample_text, sample_provenance = helpers.fetch(accession, "gsm")
        samples = [helpers.sample_fields(r) for r in helpers.parse_soft(sample_text, "SAMPLE")]
        distributions = defaultdict(Counter)
        for s in samples:
            s["series_accession"] = accession
            all_samples.append(s)
            for part in s["characteristics"].split(" | "):
                if ":" in part:
                    key, value = part.split(":", 1)
                    distributions[key.strip().lower()][value.strip()] += 1
        supplementary = [v for k, values in series.items() if k.startswith("Series_supplementary_file") for v in values]
        item = {"accession": accession, "selection_rationale": reason, "series": series,
                "series_provenance": series_provenance, "sample_provenance": sample_provenance,
                "listed_samples": len(series.get("Series_sample_id", [])), "retrieved_samples": len(samples),
                "unique_gsms": len({s["sample_accession"] for s in samples}),
                "organisms": dict(Counter(s["organism"] for s in samples)),
                "sources": dict(Counter(s["source"] for s in samples)),
                "characteristic_distributions": {k: dict(v) for k, v in distributions.items()},
                "supplementary_files": supplementary, "examples": samples[:2]}
        audit.append(item)
        (OUT / "expanded_candidate_metadata_audit.json").write_text(json.dumps(audit, ensure_ascii=False, indent=2), encoding="utf-8")
        print(json.dumps({"accession": accession, "samples": len(samples), "sources": item["sources"],
                          "keys": list(distributions), "files": supplementary,
                          "errors": [p["error"] for p in [series_provenance, sample_provenance] if "error" in p]}), flush=True)
    with (OUT / "expanded_candidate_samples.csv").open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(all_samples[0]))
        writer.writeheader()
        writer.writerows(all_samples)
    by_gsm = defaultdict(set)
    for s in all_samples:
        by_gsm[s["sample_accession"]].add(s["series_accession"])
    overlaps = [{"gsm": g, "series": sorted(v)} for g, v in by_gsm.items() if len(v) > 1]
    (OUT / "expanded_shared_gsm_audit.json").write_text(json.dumps(overlaps, indent=2), encoding="utf-8")
    print(json.dumps({"series": len(audit), "sample_records": len(all_samples), "cross_series_shared_gsms": len(overlaps)}), flush=True)


if __name__ == "__main__":
    main()
