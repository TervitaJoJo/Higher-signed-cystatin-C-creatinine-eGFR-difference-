"""Inspect cached inputs and retrieve missing official platform annotations."""
from __future__ import annotations

import gzip
import hashlib
import importlib.util
import json
import sys
import tarfile
import time
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "transcriptome_expansion_20261006"
RAW = OUT / "source_records"


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    spec = importlib.util.spec_from_file_location("download_helpers", Path(__file__).with_name("39_transcriptome_download_and_clinical_audit.py"))
    net = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(net)
    spec = importlib.util.spec_from_file_location("metadata_helpers", Path(__file__).with_name("38_transcriptome_candidate_metadata_audit.py"))
    meta = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(meta)
    ledger = []
    for acc in ["GPL15988", "GPL16686", "GPL20880"]:
        url = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?" + urlencode({"acc": acc, "targ": "self", "form": "text", "view": "full"})
        target = RAW / (acc + "_full.soft")
        if not target.exists():
            time.sleep(0.6)
            with urlopen(Request(url, headers={"User-Agent": "renal-transcriptome-expansion/1.0"}), timeout=60) as response:
                body = response.read(35_000_001)
            if len(body) > 35_000_000 or not body.startswith(b"^PLATFORM"):
                raise ValueError("Invalid or oversized public platform SOFT response")
            target.write_bytes(body)
        ledger.append({"accession": acc, "url": url, "file": str(target), "bytes": target.stat().st_size,
                       "sha256": hashlib.sha256(target.read_bytes()).hexdigest()})
    for acc in ["E-MTAB-10583", "E-MEXP-3744", "E-MTAB-9715", "E-MEXP-1937", "E-MTAB-17006", "E-MTAB-7750"]:
        url = "https://www.ebi.ac.uk/biostudies/api/v1/studies/" + acc
        ledger.append({"accession": acc, **net.get(url, RAW / (acc + "_study.json"), 3_000_000)})
    url = "https://www.ebi.ac.uk/europepmc/webservices/rest/search?" + urlencode({"query": "EXT_ID:40468979 AND SRC:MED", "format": "json", "resultType": "core"})
    ledger.append({"accession": "PMID40468979", **net.get(url, RAW / "PMID40468979_core.json", 3_000_000)})
    ledger.append({"accession": "PMC12138287", **net.get("https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12138287/fullTextXML", RAW / "PMC12138287_fulltext.xml", 8_000_000)})
    url = "https://www.ebi.ac.uk/europepmc/webservices/rest/search?" + urlencode({"query": '"GSE226151"', "format": "json", "resultType": "core"})
    ledger.append({"accession": "GSE226151_papers", **net.get(url, RAW / "GSE226151_paper_search.json", 3_000_000)})
    (OUT / "expanded_annotation_followup_ledger.json").write_text(json.dumps(ledger, indent=2), encoding="utf-8")
    if "--fetch-only" in sys.argv:
        print(json.dumps(ledger, indent=2))
        return
    samples = {}
    for f in RAW.glob("GSE*_gsm.soft"):
        samples[f.stem.split("_")[0]] = meta.parse_soft(f.read_text(encoding="utf-8"), "SAMPLE")
    for folder in sorted((OUT / "downloaded_data").iterdir()):
        acc = folder.name
        print("\nDATASET", acc)
        records = samples[acc]
        print("FIELDS", json.dumps(records[0], ensure_ascii=False)[:14000])
        if acc == "GSE15072":
            print("PLATFORM_GROUPS", [(r["accession"][0], r["Sample_platform_id"], r.get("Sample_characteristics_ch1")) for r in records])
        for f in folder.iterdir():
            print("FILE", f.name)
            if f.suffix == ".tar":
                with tarfile.open(f) as archive:
                    members = archive.getmembers()
                    print("ARCHIVE", [(m.name, m.size) for m in members][:20])
                    with gzip.open(archive.extractfile(members[0]), "rt", encoding="utf-8") as stream:
                        print("FIRST_ROWS", [next(stream).rstrip() for _ in range(3)])
            elif f.suffix == ".gz":
                with gzip.open(f, "rt", encoding="utf-8") as stream:
                    if "series_matrix" not in f.name:
                        print("FIRST_ROWS", [next(stream).rstrip()[:2800] for _ in range(3)])
                    else:
                        for line in stream:
                            if line.startswith("!Sample_data_processing"):
                                print("PROCESSING", line[:2800].rstrip())
                            if line.startswith("!series_matrix_table_begin"):
                                print("TABLE_ROWS", [next(stream).rstrip()[:2000] for _ in range(3)])
                                break
    for f in sorted(RAW.glob("GPL*")):
        if f.suffix not in {".gz", ".soft"}:
            continue
        print("\nANNOTATION", f.name)
        opener = gzip.open if f.suffix == ".gz" else open
        with opener(f, "rt", encoding="utf-8") as stream:
            for line in stream:
                if line.startswith("!platform_table_begin"):
                    print("ROWS", [next(stream).rstrip()[:2000] for _ in range(3)])
                    break
    for f in RAW.glob("E-M*_study.json"):
        text = f.read_text(encoding="utf-8")
        print("BIOSTUDIES", f.name, text[:5000])
    print("PAPER", (RAW / "PMID40468979_core.json").read_text(encoding="utf-8")[:9000])


if __name__ == "__main__":
    main()
