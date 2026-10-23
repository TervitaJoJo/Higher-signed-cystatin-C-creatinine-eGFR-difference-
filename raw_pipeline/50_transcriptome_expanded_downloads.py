"""Download phenotype-first inputs; cache provenance and cross-repository search."""
from __future__ import annotations

import importlib.util
import json
import sys
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urlencode, urljoin, urlparse

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "transcriptome_expansion_20261006"
RAW = OUT / "source_records"
DATA = OUT / "downloaded_data"
COUNTS = ["GSE111006", "GSE111010", "GSE111016", "GSE226151", "GSE174106", "GSE287726"]
ARRAYS = ["GSE117525", "GSE140358", "GSE37171", "GSE15072", "GSE142153", "GSE65218"]
PLATFORMS = ["GPL570", "GPL96", "GPL6480", "GPL15988", "GPL16686", "GPL11532"]


def helpers():
    spec = importlib.util.spec_from_file_location("expanded_downloads", Path(__file__).with_name("39_transcriptome_download_and_clinical_audit.py"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def bucket(accession: str) -> str:
    return accession[:-3] + "nnn" if len(accession[3:]) > 3 else accession[:3] + "nnn"


def main() -> None:
    sys.stdout.reconfigure(encoding="utf-8")
    DATA.mkdir(parents=True, exist_ok=True)
    net = helpers()
    audit = {x["accession"]: x for x in json.loads((OUT / "expanded_candidate_metadata_audit.json").read_text(encoding="utf-8"))}
    ledger = []

    def get(url, path, purpose, accession, cap=70_000_000):
        record = {"accession": accession, "purpose": purpose,
                  "verified_at_utc": datetime.now(timezone.utc).isoformat(), **net.get(url, path, cap)}
        ledger.append(record)
        (OUT / "expanded_download_ledger.json").write_text(json.dumps(ledger, indent=2), encoding="utf-8")
        print(json.dumps({"id": accession, "file": path.name, "bytes": record.get("bytes"), "error": record.get("error")}), flush=True)
        return record

    for accession in COUNTS:
        urls = [u.replace("ftp://", "https://", 1) for u in audit[accession]["supplementary_files"]
                if "count" in u.lower() and not u.endswith(".tar")]
        if len(urls) != 1:
            raise ValueError((accession, urls))
        url = urls[0]
        get(url, DATA / accession / Path(urlparse(url).path).name, "processed_count_matrix", accession)
    for accession in ARRAYS:
        url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{bucket(accession)}/{accession}/matrix/"
        index = RAW / f"{accession}_matrix_index.html"
        record = get(url, index, "matrix_directory_index", accession, 2_000_000)
        if record["error"]:
            continue
        parser = net.Links()
        parser.feed(index.read_text(encoding="utf-8", errors="replace"))
        links = sorted({urljoin(url, h) for h in parser.hrefs if "series_matrix" in h and h.endswith(".gz")})
        if not links:
            raise ValueError(f"No processed GEO series matrix found for {accession}")
        for source in links:
            get(source, DATA / accession / Path(urlparse(source).path).name, "processed_array_series_matrix", accession)
    accession = "GSE129534"
    url = audit[accession]["supplementary_files"][0].replace("ftp://", "https://", 1)
    get(url, DATA / accession / Path(urlparse(url).path).name, "per_sample_processed_counts_archive", accession)

    for platform in PLATFORMS:
        url = f"https://ftp.ncbi.nlm.nih.gov/geo/platforms/{bucket(platform)}/{platform}/annot/{platform}.annot.gz"
        get(url, RAW / f"{platform}.annot.gz", "official_probe_annotation", platform, 30_000_000)

    cross = []
    for query in ["sarcopenia", "frailty", '"chronic kidney disease" AND muscle']:
        params = {"query": query, "collection": "ArrayExpress", "pageSize": 200, "page": 1}
        url = "https://www.ebi.ac.uk/biostudies/api/v1/search?" + urlencode(params)
        target = RAW / ("biostudies_" + ("ckd_muscle" if "kidney" in query else query) + ".json")
        record = get(url, target, "cross_repository_search", "BioStudies_ArrayExpress", 10_000_000)
        if record["error"]:
            cross.append({"query": query, **record})
            continue
        content = json.loads(target.read_text(encoding="utf-8"))
        hits = content.get("hits", [])
        cross.append({"query": query, "totalHits": content.get("totalHits"), "retrieved": len(hits),
                      "complete_for_query": content.get("totalHits") == len(hits),
                      "non_geo_accessions": [h["accession"] for h in hits if not h["accession"].startswith("E-GEOD-")],
                      "human_relevant_hits": [h for h in hits if "homo sapiens" in h.get("content", "").lower()],
                      "file": str(target), "url": url})
    (OUT / "biostudies_cross_repository_audit.json").write_text(json.dumps(cross, indent=2), encoding="utf-8")
    print(json.dumps({"completed_download_records": len(ledger), "download_errors": sum(bool(x.get("error")) for x in ledger)}), flush=True)


if __name__ == "__main__":
    main()
