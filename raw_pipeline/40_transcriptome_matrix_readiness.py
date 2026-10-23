"""Create auditable matrix-to-GSM links, QC and a supplement evidence check."""
from __future__ import annotations

import csv
import hashlib
import importlib.util
import io
import json
import re
import sys
import zipfile
import xml.etree.ElementTree as ET
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
import pymupdf

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "transcriptome_feasibility"
PREP = OUT / "analysis_ready"


def load_download_module():
    spec = importlib.util.spec_from_file_location("download_audit", Path(__file__).with_name("39_transcriptome_download_and_clinical_audit.py"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def characteristics(text: str) -> dict:
    result = {}
    for part in text.split(" | "):
        if ":" in part:
            key, value = part.split(":", 1)
            result[key.strip().lower()] = value.strip()
    return result


def sample_column(accession: str, row: dict) -> str | None:
    if accession == "GSE236927":
        return row["sample_accession"]
    if accession in {"GSE157712", "GSE309619"}:
        match = re.search(r"\[([^\]]+)\]$", row["title"])
        return match.group(1) if match else None
    if accession == "GSE287118":
        return row["title"]
    if accession == "GSE144304":
        return row["title"].split("_", 1)[0]
    if accession == "GSE339456":
        prefixes = {"EarlyStable": "ES", "EarlyProgressive": "EP", "LateStable": "LS", "LateProgressive": "LP"}
        match = re.search(r",\s*(\w+),\s*nr(\d+)", row["title"])
        return prefixes[match.group(1)] + "_" + match.group(2) if match else None
    return None


def supplement_audit(download_module) -> list[dict]:
    reports = []
    for pmcid in ["PMC10877806", "PMC8783691"]:
        url = f"https://www.ebi.ac.uk/europepmc/webservices/rest/{pmcid}/supplementaryFiles"
        path = OUT / "source_records" / (pmcid + "_supplementaryFiles.zip")
        if pmcid == "PMC8783691":
            url = "https://insight.jci.org/articles/view/150112/sd/pdf/render/1"
            path = OUT / "source_records" / "PMC8783691_jciinsight-6-150112-s076.pdf"
        record = {"pmcid": pmcid, **download_module.get(url, path, 20_000_000)}
        if pmcid == "PMC8783691":
            record["prior_zip_url"] = f"https://www.ebi.ac.uk/europepmc/webservices/rest/{pmcid}/supplementaryFiles"
            record["prior_zip_error"] = "Download exceeds 20000000 bytes;used individually linked publisher supplement"
        if record["error"]:
            reports.append(record)
            continue
        record["members"] = []
        if path.suffix == ".pdf":
            payloads = [(path.name, path.read_bytes())]
        else:
            with zipfile.ZipFile(path) as archive:
                payloads = [(info.filename, archive.read(info)) for info in archive.infolist()
                            if info.file_size <= 25_000_000 and info.filename.lower().endswith((".pdf", ".docx"))]
        for name, payload in payloads:
                member = {"name": name, "bytes": len(payload)}
                if name.lower().endswith(".pdf"):
                    with pymupdf.open(stream=payload, filetype="pdf") as doc:
                        texts = [page.get_text() for page in doc]
                    member["pages"] = len(texts)
                    snippets = []
                    for i, text in enumerate(texts, 1):
                        for match in re.finditer(r"cystatin|creatinine|eGFR|estimated glomerular|GSE\d+|subject|individual|RSEM|salmon", text, re.I):
                            snippets.append({"page": i, "text": text[max(0, match.start()-100):match.end()+280].replace("\n", " ")})
                    member["clinical_or_processing_excerpts"] = snippets[:50]
                    member["text"] = "\n".join(texts)
                else:
                    with zipfile.ZipFile(io.BytesIO(payload)) as docx:
                        document = ET.fromstring(docx.read("word/document.xml"))
                    ns = {"w": "http://schemas.openxmlformats.org/wordprocessingml/2006/main"}
                    paragraphs = ["".join(p.itertext()) for p in document.findall(".//w:p", ns)]
                    member["text"] = "\n".join(paragraphs)
                # Preserve extraction for inspection; do not execute embedded code.
                stem = Path(name).name
                target = OUT / "source_records" / (pmcid + "_" + stem + ".txt")
                target.write_text(member.pop("text"), encoding="utf-8")
                member["extracted_text_file"] = str(target)
                record["members"].append(member)
        reports.append(record)
        print(json.dumps({"pmcid": pmcid, "supplement_members": len(record["members"]), "error": record["error"]}), flush=True)
    return reports


def main() -> None:
    sys.stdout.reconfigure(encoding="utf-8")
    PREP.mkdir(parents=True, exist_ok=True)
    if "--supplements-only" in sys.argv:
        supplements = supplement_audit(load_download_module())
        (OUT / "clinical_supplement_audit.json").write_text(json.dumps(supplements, ensure_ascii=False, indent=2), encoding="utf-8")
        return
    with (OUT / "candidate_sample_metadata.csv").open(encoding="utf-8-sig", newline="") as handle:
        samples = list(csv.DictReader(handle))
    manifest = json.loads((OUT / "download_manifest_and_matrix_audit.json").read_text(encoding="utf-8"))
    results = []
    for entry in manifest:
        accession = entry["accession"]
        if entry.get("error") or ".R.gz" in entry["file"]:
            continue
        path = Path(entry["file"])
        sep = ";" if accession == "GSE339456" else "," if accession == "GSE287118" else "\t"
        frame = pd.read_csv(path, sep=sep, comment="!", index_col=0, low_memory=False)
        metadata = [s for s in samples if s["series_accession"] == accession]
        mapping = []
        missing = []
        for sample in metadata:
            token = sample_column(accession, sample)
            matches = [c for c in frame.columns if c == token or (accession == "GSE144304" and c.startswith(str(token) + "_"))]
            if len(matches) != 1:
                missing.append({"gsm": sample["sample_accession"], "token": token, "matches": matches})
                continue
            clinical = characteristics(sample["characteristics"])
            mapping.append({"sample_accession": sample["sample_accession"], "matrix_column": matches[0], "sample_title": sample["title"], **clinical})
        columns = [m["matrix_column"] for m in mapping]
        unique_mapping = len(columns) == len(set(columns))
        numeric = frame.loc[:, columns].apply(pd.to_numeric, errors="coerce")
        values = numeric.to_numpy(dtype=float)
        finite = np.isfinite(values)
        expected_count_type = accession in {"GSE157712", "GSE144304", "GSE339456"}
        donor_ids = [m["subject id"] for m in mapping if "subject id" in m]
        repeated = {key: value for key, value in Counter(donor_ids).items() if value > 1}
        result = {
            "accession": accession, "source_matrix": str(path), "features": len(frame),
            "source_columns": len(frame.columns), "matched_sample_columns": len(columns),
            "metadata_sample_records": len(metadata), "missing_sample_links": missing,
            "mapping_is_one_to_one": unique_mapping, "duplicate_feature_ids": int(frame.index.duplicated().sum()),
            "nonfinite_expression_entries": int((~finite).sum()),
            "negative_expression_entries": int((values[finite] < 0).sum()),
            "fractional_entry_fraction": float(np.mean(np.abs(values[finite] - np.rint(values[finite])) > 1e-7)),
            "count_type_from_source": expected_count_type,
            "unique_subject_ids_if_available": len(set(donor_ids)) if donor_ids else None,
            "repeated_subject_ids": repeated,
            "unmatched_source_columns": [str(c) for c in frame.columns if c not in columns],
            "ready_for_genomewide_model": bool(len(frame) > 1000 and not missing and unique_mapping and finite.all()),
        }
        if accession == "GSE157712":
            result["analysis_constraint"] = "RSEM expected counts include fractional values; do not silently round. Primary CKD vs control: 7+7 independent donors. Dialysis sensitivity must account for repeated subject IDs 10,14,15."
        elif accession == "GSE309619":
            result["analysis_constraint"] = "Only four features in public FPKM file. Twenty GSMs are five donors x two visits x two biopsies/replicates. One additional unmatched expression column. Genomewide analysis requires obtaining full counts or rebuilding from raw reads."
            result["derived_donors"] = len(set(re.findall(r"Patient GHD(\d+)", " | ".join(s["title"] for s in metadata))))
        elif accession == "GSE144304":
            result["analysis_constraint"] = "Primary functional comparison is older adults only (24 Frail,30 Fit); 26 Young samples are secondary aging comparators. Use original group, not derived logitgroup."
        elif accession == "GSE236927":
            result["analysis_constraint"] = "95 women; FI available in 88. No individual serum cystatin C/creatinine in GEO characteristics. Published monocyte/frailty analysis already exists."
        elif accession == "GSE287118":
            result["analysis_constraint"] = "Only TPM in public processed file. Four disease groups; PAD confounding must be modeled. No public individual renal/function covariates in GEO characteristics."
        else:
            result["analysis_constraint"] = "35 biopsy samples; group obtainable from title. No individual eGFR trajectories/covariates in GEO characteristics. Supports renal progression, not ADL recovery."
        fields = list(dict.fromkeys(k for item in mapping for k in item))
        with (PREP / (accession + "_sample_sheet.csv")).open("w", encoding="utf-8-sig", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields)
            writer.writeheader(); writer.writerows(mapping)
        if result["ready_for_genomewide_model"]:
            numeric.to_csv(PREP / (accession + "_expression.tsv.gz"), sep="\t", compression="gzip")
        results.append(result)
        print(json.dumps({k:result[k] for k in ["accession", "features", "matched_sample_columns", "missing_sample_links", "unique_subject_ids_if_available", "ready_for_genomewide_model"]}), flush=True)
    (OUT / "matrix_readiness_audit.json").write_text(json.dumps(results, ensure_ascii=False, indent=2), encoding="utf-8")
    module = load_download_module()
    supplements = supplement_audit(module)
    (OUT / "clinical_supplement_audit.json").write_text(json.dumps(supplements, ensure_ascii=False, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
