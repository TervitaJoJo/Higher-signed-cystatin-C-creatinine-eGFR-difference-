"""Build audited matrices and sample sheets without phenotype guesses."""
from __future__ import annotations

import gzip
import hashlib
import importlib.util
import io
import json
import re
import sys
import tarfile
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "transcriptome_expansion_20261006"
RAW = OUT / "source_records"
READY = OUT / "analysis_ready"


def embedded_table(path, start, end):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as h:
        active = False
        lines = []
        for line in h:
            if line.startswith(start):
                active = True
            elif line.startswith(end):
                break
            elif active:
                lines.append(line)
    if not lines:
        raise ValueError(f"No embedded table: {path}")
    return pd.read_csv(io.StringIO("".join(lines)), sep="\t", dtype={"ID": str, "ID_REF": str})


def metadata(record):
    fields = {}
    for part in record.get("Sample_characteristics_ch1", []):
        if ":" in part:
            k, v = part.split(":", 1)
            k = k.strip().lower()
            if k in fields:
                raise ValueError((record["accession"], k))
            fields[k] = v.strip()
    row = {"sample_accession": record["accession"][0], "title": record["Sample_title"][0],
           "platform": record["Sample_platform_id"][0],
           "source": record["Sample_source_name_ch1"][0],
           "original_characteristics": " | ".join(record.get("Sample_characteristics_ch1", []))}
    row.update(fields)
    row["sex"] = fields.get("sex", "").lower()
    row["age"] = next((fields[k] for k in ["age", "age (yr)", "age (yrs)", "age (years)", "age_tx"] if k in fields), "")
    row["bmi"] = next((fields[k] for k in ["bmi", "bmi (kg/m2)", "body mass index"] if k in fields), "")
    return row


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    READY.mkdir(parents=True, exist_ok=True)
    spec = importlib.util.spec_from_file_location("metadata_helpers", Path(__file__).with_name("38_transcriptome_candidate_metadata_audit.py"))
    meta = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(meta)
    audit = {}
    contracts = []
    for folder in sorted((OUT / "downloaded_data").iterdir()):
        acc = folder.name
        if acc == "GSE287726":
            audit[acc] = {"excluded_from_models": "No public individual frailty/clinical mapping; muscle, not blood"}
            continue
        records = meta.parse_soft((RAW / f"{acc}_gsm.soft").read_text(encoding="utf-8"), "SAMPLE")
        rows = [metadata(r) for r in records]
        if acc.startswith("GSE1110"):
            for row in rows:
                row["matrix_key"] = re.fullmatch(r"(Sample \d+) \[[a-z]+\]", row["title"])[1]
                yes = row["sarcopenia status"] == "yes"
                partial = row.get("low muscle mass") == "yes" or row.get("low muscle strength and/or low physical performance") == "yes"
                row["group"] = "Sarcopenia" if yes else ("Partial" if partial else "Healthy")
        elif acc == "GSE226151":
            for row in rows:
                label, num = row["title"].split(", ")[1].split(" ")
                row["matrix_key"] = label + "_" + str(int(num))
                row["group"] = {"healthy aged": "Healthy", "pre-sarcopenia": "Presarcopenia", "sarcopenia": "Sarcopenia"}[row["disease state"]]
        elif acc == "GSE174106":
            for row in rows:
                row["matrix_key"] = row["title"]
                row["group"] = {"sarcopenic": "ALMloss", "stable high": "Stablehigh"}.get(row["appendicular lean mass (alm)"], "Other")
        elif acc == "GSE129534":
            for row in rows:
                row["matrix_key"] = row["sample_accession"]
                row["group"] = {"Frail": "Frail", "Non-Frail": "Nonfrail"}[row["frailty status"]]
        elif acc == "GSE117525":
            for row in rows:
                row["matrix_key"] = row["sample_accession"]
                row["group"] = {"Frail (baseline)": "Frail", "Healthy older (baseline)": "Healthy"}.get(row["group"], "Other")
                row["subject_id"] = row["subjectid"]
        elif acc == "GSE65218":
            for row in rows:
                row["matrix_key"] = row["sample_accession"]
                row["group"] = {"frail": "Frail", "non-frail": "Nonfrail", "pre-frail": "Prefrail"}[row["frailty index"]]
                row["subject_id"] = row["title"].split(" ")[-1]
        elif acc == "GSE140358":
            for row in rows:
                row["matrix_key"] = row["sample_accession"]
                row["group"] = row["condition"]
        elif acc == "GSE37171":
            for row in rows:
                row["matrix_key"] = row["sample_accession"]
                row["group"] = {"uremia": "Uremia", "healthy control": "Healthy"}[row["disease state"]]
                row["study_set"] = re.search(r"\[([^]]+)\]", row["title"])[1]
                row["subject_id"] = row["title"].split(" [")[0]
        elif acc == "GSE15072":
            for row in rows:
                row["matrix_key"] = row["sample_accession"]
                row["group"] = "Healthy" if "healthy subject" in row else ("CKD" if row["treatment"] == "conservative" else "HD")
                row["subject_id"] = row.get("patient", row.get("healthy subject", ""))
        elif acc == "GSE142153":
            for row in rows:
                row["matrix_key"] = row["sample_accession"]
                row["group"] = {"healthy control": "Healthy", "diabetic nephropathy": "DN", "ESRD": "ESRD"}[row["diagnosis"]]
        else:
            raise ValueError(acc)
        sheets = pd.DataFrame(rows).set_index("matrix_key", drop=False)
        if not sheets.index.is_unique:
            raise ValueError(f"Ambiguous sample mapping: {acc}")
        for f in sorted(folder.iterdir()):
            if f.suffix == ".tar":
                columns = []
                with tarfile.open(f) as archive:
                    for member in archive.getmembers():
                        if not member.isfile() or not re.fullmatch(r"GSM\d+_[A-Za-z0-9_]+\.count\.txt\.gz", member.name):
                            raise ValueError(f"Unexpected archive entry: {member.name}")
                        with gzip.open(archive.extractfile(member), "rt", encoding="utf-8") as stream:
                            tab = pd.read_csv(stream, sep="\t", index_col=0)
                        assert tab.shape[1] == 1
                        gsm = member.name.split("_")[0]
                        assert tab.columns[0] == "(count) " + sheets.loc[gsm, "title"]
                        tab.columns = [gsm]
                        columns.append(tab)
                assert all(c.index.equals(columns[0].index) for c in columns)
                x = pd.concat(columns, axis=1)
            elif "series_matrix" in f.name:
                x = embedded_table(f, "!series_matrix_table_begin", "!series_matrix_table_end").set_index("ID_REF")
            else:
                x = pd.read_csv(f, sep="," if f.name.endswith(".csv.gz") else "\t", index_col=0)
            x.index = x.index.astype(str)
            x.columns = x.columns.astype(str)
            keys = [c.split("|")[0] if acc == "GSE174106" else c for c in x.columns]
            assert len(set(keys)) == len(keys) and set(keys) <= set(sheets.index)
            m = sheets.loc[keys].reset_index(drop=True)
            m.insert(0, "matrix_column", x.columns)
            assert not m["sample_accession"].duplicated().any()
            assert x.index.is_unique and np.isfinite(x.to_numpy()).all()
            platforms = m["platform"].unique()
            assert len(platforms) == 1
            platform = platforms[0]
            name = acc + ("_" + platform if acc == "GSE15072" else "")
            pairs = []
            for i, col in enumerate(x.columns):
                for j in range(i + 1, x.shape[1]):
                    if np.array_equal(x.iloc[:, i].values, x.iloc[:, j].values):
                        pairs.append([col, x.columns[j]])
            counts = "series_matrix" not in f.name
            if counts:
                assert (x.to_numpy() >= 0).all()
                assay = "counts"
                keytype = "ENSEMBL" if acc.startswith("GSE1110") else ("SYMBOL_ENTREZ" if acc == "GSE174106" else "SYMBOL")
            else:
                assay = "array_MAS5_linear" if acc == "GSE15072" else "array_deposited_log"
                keytype = "ENTREZ_SUFFIX" if acc == "GSE117525" else "PROBE"
            x.index.name = "feature_id"
            x.to_csv(READY / (name + "_expression.tsv.gz"), sep="\t", compression="gzip")
            m.to_csv(READY / (name + "_sample_sheet.csv"), index=False, encoding="utf-8-sig")
            entry = {"accession": acc, "name": name, "platform": platform, "assay": assay, "keytype": keytype,
                     "features": x.shape[0], "samples": x.shape[1], "value_quantiles": np.quantile(x.to_numpy(), [0, .01, .5, .99, 1]).tolist(),
                     "groups": m["group"].value_counts().to_dict(), "exact_duplicate_columns": pairs,
                     "input_sha256": hashlib.sha256(f.read_bytes()).hexdigest(),
                     "provisional_only": acc == "GSE226151"}
            if acc == "GSE226151":
                entry["warning"] = "Human sample metadata and symbols but all sample processing says mm10/mouse and RPKM despite filename READ_COUNT; exclude from qualified PCR evidence until resolved"
            audit[name] = entry
            contracts.append(entry)
            print(json.dumps(entry), flush=True)
    for platform in ["GPL570", "GPL96", "GPL6480", "GPL15988", "GPL16686", "GPL20880"]:
        f = RAW / (platform + ("_full.soft" if platform in {"GPL15988", "GPL16686", "GPL20880"} else ".annot.gz"))
        table = embedded_table(f, "!platform_table_begin", "!platform_table_end").fillna("")
        assert not table["ID"].duplicated().any()
        ann = pd.DataFrame({"feature_id": table["ID"]})
        ann["source_symbol"] = table.get("Gene symbol", table.get("GENE SYMBOL", pd.Series("", index=table.index)))
        ann["source_entrez"] = table.get("Gene ID", pd.Series("", index=table.index))
        ann["source_refseq"] = table.get("GB_ACC", pd.Series("", index=table.index))
        ann.to_csv(READY / (platform + "_annotation.csv"), index=False, encoding="utf-8-sig")
        print("ANNOTATION", platform, table.shape, list(table.columns), flush=True)
    (OUT / "expanded_input_audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    (READY / "analysis_contracts.json").write_text(json.dumps(contracts, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
