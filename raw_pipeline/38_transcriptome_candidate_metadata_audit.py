"""Retrieve series and sample SOFT records for a manually scoped shortlist."""
from __future__ import annotations

import csv
import json
import re
import sys
import time
from collections import Counter
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "transcriptome_feasibility"
RAW = OUT / "source_records"
CANDIDATES = {
    "GSE236927": "older PBMC; cystatin C and frailty study",
    "GSE68526": "HRS blood RNA pilot; verify public linkage and renal phenotypes",
    "GSE157712": "CKD muscle; kidney function and physical function study",
    "GSE287118": "CKD and PAD calf muscle; check factorial confounding and donors",
    "GSE309619": "dialysis muscle before/after exercise; verify pairing",
    "GSE144304": "fit and pre-frail muscle; independent functional comparator",
    "GSE111017": "sarcopenia muscle superseries; verify constituent studies",
    "GSE142025": "diabetic kidney biopsy RNA-seq; clinical eGFR support",
    "GSE104954": "ERCB tubulointerstitium; explicitly check previously reused samples",
    "GSE66494": "independent CKD biopsy microarray; clinical phenotype support",
    "GSE339456": "hypertensive nephropathy progression; verify public follow-up data",
    "GSE183276": "human kidney cell-state localization; not primary functional validation",
}


def fetch(accession: str, target: str) -> tuple[str, dict]:
    if not re.fullmatch(r"GSE\d+", accession) or target not in {"self", "gsm"}:
        raise ValueError((accession, target))
    path = RAW / f"{accession}_{target}.soft"
    params = {"acc": accession, "targ": target, "form": "text", "view": "brief"}
    url = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?" + urlencode(params)
    if path.exists():
        return path.read_text(encoding="utf-8"), {"url": url, "cached": True}
    error = None
    for attempt in range(2):
        try:
            time.sleep(0.6)
            req = Request(url, headers={"User-Agent": "renal-transcriptome-feasibility/1.0"})
            with urlopen(req, timeout=50) as response:
                content = response.read(12_000_001)
            if len(content) > 12_000_000:
                raise RuntimeError("SOFT exceeds the bounded metadata budget")
            text = content.decode("utf-8", errors="replace")
            if f"^SERIES = {accession}" not in text and "^SAMPLE = GSM" not in text:
                raise RuntimeError("Response is not the requested GEO SOFT metadata")
            path.write_text(text, encoding="utf-8")
            return text, {"url": url, "cached": False, "bytes": len(content)}
        except Exception as exc:
            error = str(exc)
            if attempt == 0:
                time.sleep(3)
    return "", {"url": url, "error": error}


def parse_soft(text: str, kind: str) -> list[dict[str, list[str]]]:
    records = []
    active = None
    for line in text.splitlines():
        if line.startswith("^"):
            if active is not None:
                records.append(active)
            active = {"accession": [line.split(" = ", 1)[1]]} if line.startswith(f"^{kind} = ") else None
        elif active is not None and line.startswith("!") and " = " in line:
            key, value = line[1:].split(" = ", 1)
            active.setdefault(key, []).append(value)
    if active is not None:
        records.append(active)
    return records


def sample_fields(record: dict) -> dict:
    characteristics = record.get("Sample_characteristics_ch1", [])
    return {
        "sample_accession": record["accession"][0],
        "title": "; ".join(record.get("Sample_title", [])),
        "source": "; ".join(record.get("Sample_source_name_ch1", [])),
        "organism": "; ".join(record.get("Sample_organism_ch1", [])),
        "platform": "; ".join(record.get("Sample_platform_id", [])),
        "characteristics": " | ".join(characteristics),
        "relations": " | ".join(record.get("Sample_relation", [])),
        "supplementary_files": " | ".join(v for k, vv in record.items() if k.startswith("Sample_supplementary_file") for v in vv),
    }


def main() -> None:
    sys.stdout.reconfigure(encoding="utf-8")
    RAW.mkdir(parents=True, exist_ok=True)
    audits = []
    all_samples = []
    for accession, rationale in CANDIDATES.items():
        series_text, series_provenance = fetch(accession, "self")
        series_records = parse_soft(series_text, "SERIES")
        series = series_records[0] if series_records else {}
        sample_text, sample_provenance = fetch(accession, "gsm")
        samples = [sample_fields(x) for x in parse_soft(sample_text, "SAMPLE")]
        characteristic_counts: Counter = Counter()
        for sample in samples:
            sample["series_accession"] = accession
            all_samples.append(sample)
            for part in sample["characteristics"].split(" | "):
                if ":" in part:
                    characteristic_counts[part.split(":", 1)[0].strip().lower()] += 1
        supplementary = [v for k, vv in series.items() if k.startswith("Series_supplementary_file") for v in vv]
        audit = {
            "accession": accession, "rationale": rationale,
            "series_record": series, "series_provenance": series_provenance,
            "sample_provenance": sample_provenance, "retrieved_sample_records": len(samples),
            "series_listed_samples": len(series.get("Series_sample_id", [])),
            "characteristic_key_counts": dict(characteristic_counts),
            "sample_sources": dict(Counter(s["source"] for s in samples)),
            "supplementary_files": supplementary,
            "example_samples": samples[:3],
        }
        audits.append(audit)
        print(json.dumps({"accession": accession, "sample_records": len(samples),
                          "characteristic_keys": list(characteristic_counts),
                          "series_supplementary_files": supplementary,
                          "errors": [x.get("error") for x in (series_provenance, sample_provenance) if x.get("error")]}), flush=True)
    if all_samples:
        with (OUT / "candidate_sample_metadata.csv").open("w", encoding="utf-8-sig", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(all_samples[0]))
            writer.writeheader()
            writer.writerows(all_samples)
    (OUT / "candidate_metadata_audit.json").write_text(json.dumps(audits, ensure_ascii=False, indent=2), encoding="utf-8")
    overlap = []
    series_by_gsm = {}
    for sample in all_samples:
        series_by_gsm.setdefault(sample["sample_accession"], set()).add(sample["series_accession"])
    for gsm, series_set in series_by_gsm.items():
        if len(series_set) > 1:
            overlap.append({"gsm": gsm, "series": sorted(series_set)})
    (OUT / "shared_gsm_audit.json").write_text(json.dumps(overlap, indent=2), encoding="utf-8")
    print(json.dumps({"candidate_series": len(audits), "sample_records": len(all_samples), "shared_gsms": len(overlap)}), flush=True)


if __name__ == "__main__":
    main()
