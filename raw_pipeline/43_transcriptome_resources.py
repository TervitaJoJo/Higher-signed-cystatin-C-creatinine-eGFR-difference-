"""Cache small public pathway libraries and official Illumina annotations."""
from __future__ import annotations

import csv
import gzip
import hashlib
import io
import json
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

import requests

OUT = ROOT / "transcriptome_analysis" / "resources"
LIBRARIES = ("MSigDB_Hallmark_2020", "Reactome_2022")


def retrieve(url: str, target: Path, cap: int = 30_000_000) -> dict:
    if not target.exists():
        last = None
        for attempt in range(3):
            try:
                with requests.get(url, timeout=(20, 120), stream=True) as response:
                    response.raise_for_status()
                    payload = bytearray()
                    for chunk in response.iter_content(131072):
                        payload.extend(chunk)
                        if len(payload) > cap:
                            raise ValueError(f"Resource exceeds {cap} bytes: {url}")
                    target.write_bytes(payload)
                break
            except requests.RequestException as exc:
                last = exc
                if attempt == 2:
                    raise
                time.sleep(2 * (attempt + 1))
        if not target.exists():
            raise RuntimeError(str(last))
    payload = target.read_bytes()
    return {"url": url, "file": str(target), "bytes": len(payload),
            "sha256": hashlib.sha256(payload).hexdigest(),
            "verified_at_utc": datetime.now(timezone.utc).isoformat()}


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    ledger = []
    stats_url = "https://maayanlab.cloud/Enrichr/datasetStatistics"
    ledger.append(retrieve(stats_url, OUT / "enrichr_datasetStatistics.json"))
    names = {x["libraryName"] for x in json.loads((OUT / "enrichr_datasetStatistics.json").read_text())["statistics"]}
    for name in LIBRARIES:
        if name not in names:
            raise ValueError(f"Library not in current Enrichr catalog: {name}")
        url = f"https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName={name}"
        record = retrieve(url, OUT / f"{name}.gmt")
        lines = (OUT / f"{name}.gmt").read_text(encoding="utf-8").splitlines()
        if not lines or any(len(line.split("\t")) < 3 for line in lines):
            raise ValueError(f"Malformed GMT: {name}")
        record["library"] = name
        record["gene_sets"] = len(lines)
        ledger.append(record)
        print(json.dumps({"library": name, "gene_sets": len(lines)}), flush=True)

    ledger.append({"url": "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL10nnn/GPL10558/soft/GPL10558_family.soft.gz",
                   "status": "prior_attempt_stopped_at_30000000_byte_cap",
                   "reason": "Platform family is larger than required; retrieve platform self only."})
    ledger.append({"url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL10558&targ=self&form=text&view=full",
                   "status": "prior_attempt_stopped_at_30000000_byte_cap",
                   "reason": "Full platform includes sequence columns; retrieve the annotation-only file."})
    url = "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL10nnn/GPL10558/annot/GPL10558.annot.gz"
    path = OUT / "GPL10558.annot.gz"
    record = retrieve(url, path)
    text = gzip.decompress(path.read_bytes()).decode("utf-8", errors="replace")
    table = text.split("!platform_table_begin", 1)[1].split("!platform_table_end", 1)[0].strip()
    reader = csv.DictReader(io.StringIO(table), delimiter="\t")
    rows = list(reader)
    required = {"ID", "Gene symbol", "Gene ID"}
    if not required.issubset(reader.fieldnames or []):
        raise ValueError(f"Unexpected GPL columns: {reader.fieldnames}")
    with (OUT / "GPL10558_probe_annotation.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["feature_id", "source_symbol", "entrez_id"])
        writer.writeheader()
        for row in rows:
            writer.writerow({"feature_id": row["ID"], "source_symbol": row["Gene symbol"],
                             "entrez_id": row["Gene ID"]})
    record["platform_features"] = len(rows)
    ledger.append(record)
    (OUT / "retrieval_ledger.json").write_text(json.dumps(ledger, indent=2), encoding="utf-8")
    print(json.dumps({"platform": "GPL10558", "features": len(rows)}), flush=True)


if __name__ == "__main__":
    main()
