"""Download bounded processed data and audit matching clinical metadata."""
from __future__ import annotations

import csv
import gzip
import hashlib
import json
import re
import sys
import time
import xml.etree.ElementTree as ET
from collections import Counter
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import urlencode, urljoin, urlparse
from urllib.request import Request, urlopen

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "transcriptome_feasibility"
DATA = OUT / "downloaded_data"
RAW = OUT / "source_records"
LIMIT = 70_000_000
FILES = {
    "GSE236927": [
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE236nnn/GSE236927/matrix/GSE236927_series_matrix.txt.gz",
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE236nnn/GSE236927/suppl/GSE236927_run_GCF-2016-413.R.gz",
    ],
    "GSE157712": ["https://ftp.ncbi.nlm.nih.gov/geo/series/GSE157nnn/GSE157712/suppl/GSE157712_Muscle_CKD_raw_counts.txt.gz"],
    "GSE287118": ["https://ftp.ncbi.nlm.nih.gov/geo/series/GSE287nnn/GSE287118/suppl/GSE287118_TPM_values.csv.gz"],
    "GSE144304": ["https://ftp.ncbi.nlm.nih.gov/geo/series/GSE144nnn/GSE144304/suppl/GSE144304_JoinedCounts.txt.gz"],
    "GSE309619": ["https://ftp.ncbi.nlm.nih.gov/geo/series/GSE309nnn/GSE309619/suppl/GSE309619_AR_FPKM.xls.txt.gz"],
    "GSE339456": ["https://ftp.ncbi.nlm.nih.gov/geo/series/GSE339nnn/GSE339456/suppl/GSE339456_raw_gene_count.csv.gz"],
}
PMIDS = ["38373890", "34784301", "39963788", "36820956", "37430322", "37191430", "41831133", "26246388"]


def get(url: str, target: Path, limit: int = LIMIT) -> dict:
    parsed = urlparse(url)
    if parsed.scheme != "https" or parsed.hostname not in {"ftp.ncbi.nlm.nih.gov", "www.ebi.ac.uk", "insight.jci.org"}:
        raise ValueError("Only authoritative public download hosts are permitted")
    error = None
    if not target.exists():
        for attempt in range(2):
            try:
                time.sleep(0.6)
                req = Request(url, headers={"User-Agent": "renal-transcriptome-feasibility/1.0"})
                with urlopen(req, timeout=60) as response:
                    length = int(response.headers.get("Content-Length", "0"))
                    if length > limit:
                        raise RuntimeError(f"Download exceeds {limit} bytes ({length})")
                    body = response.read(limit + 1)
                if len(body) > limit:
                    raise RuntimeError(f"Download exceeds {limit} bytes")
                target.parent.mkdir(parents=True, exist_ok=True)
                target.write_bytes(body)
                break
            except Exception as exc:
                error = str(exc)
                if attempt == 0:
                    time.sleep(3)
    if not target.exists():
        return {"url": url, "file": str(target), "error": error}
    return {"url": url, "file": str(target), "bytes": target.stat().st_size,
            "sha256": hashlib.sha256(target.read_bytes()).hexdigest(), "error": None}


def matrix_preview(path: Path) -> dict:
    opener = gzip.open if path.suffix == ".gz" else open
    header = None
    rows = 0
    first_rows = []
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("!") or not line.strip():
                continue
            if header is None:
                delimiter = csv.Sniffer().sniff(line, delimiters="\t,;").delimiter
                header = next(csv.reader([line], delimiter=delimiter))
                continue
            values = next(csv.reader([line], delimiter=delimiter))
            rows += 1
            if len(first_rows) < 2:
                first_rows.append(values[:6])
    return {"header": header, "data_rows": rows, "example_values": first_rows}


class Links(HTMLParser):
    def __init__(self):
        super().__init__()
        self.hrefs = []

    def handle_starttag(self, tag, attrs):
        if tag.lower() == "a":
            self.hrefs.extend(v for k, v in attrs if k.lower() == "href" and v)


def summarize_samples() -> list[dict]:
    with (OUT / "candidate_sample_metadata.csv").open(encoding="utf-8-sig", newline="") as handle:
        samples = list(csv.DictReader(handle))
    summaries = []
    for accession in sorted({s["series_accession"] for s in samples}):
        selected = [s for s in samples if s["series_accession"] == accession]
        by_key = {}
        for sample in selected:
            for part in sample["characteristics"].split(" | "):
                if ":" in part:
                    key, value = part.split(":", 1)
                    by_key.setdefault(key.strip().lower(), []).append(value.strip())
        cats = {key: dict(Counter(values)) for key, values in by_key.items()
                if any(term in key for term in ("group", "state", "status", "condition", "sex", "gender", "population", "time point", "diagnosis", "study set"))}
        numeric = {}
        for key, values in by_key.items():
            if any(term in key for term in ("age", "bmi", "grip", "sppb", "walktime", "fi score", "egfr", "cystatin", "creatinine")):
                parsed = []
                for value in values:
                    try:
                        number = float(value)
                        if number == number and abs(number) != float("inf"):
                            parsed.append(number)
                    except ValueError:
                        pass
                numeric[key] = {"n_numeric": len(parsed), "min": min(parsed) if parsed else None, "max": max(parsed) if parsed else None}
        summaries.append({"accession": accession, "sample_records": len(selected), "categorical_distributions": cats,
                          "numeric_phenotypes": numeric, "characteristic_keys": list(by_key),
                          "unique_subject_ids": len(set(by_key.get("subject id", []))) if "subject id" in by_key else None,
                          "unique_patient_ids": len(set(by_key.get("patient", []))) if "patient" in by_key else None})
    return summaries


def paper_audit() -> list[dict]:
    query = " OR ".join("EXT_ID:" + x for x in PMIDS)
    url = "https://www.ebi.ac.uk/europepmc/webservices/rest/search?" + urlencode({"query": query, "format": "json", "resultType": "core", "pageSize": 25})
    result = get(url, RAW / "candidate_papers_epmc.json", 5_000_000)
    if result["error"]:
        return [result]
    records = json.loads((RAW / "candidate_papers_epmc.json").read_text(encoding="utf-8"))
    papers = []
    for record in records.get("resultList", {}).get("result", []):
        paper = {"pmid": record.get("id"), "pmcid": record.get("pmcid"), "title": record.get("title"),
                 "doi": record.get("doi"), "year": record.get("pubYear"), "is_open_access": record.get("isOpenAccess"),
                 "source_query": url, "supplementary_materials": [], "data_availability": [], "clinical_excerpts": []}
        pmcid = paper["pmcid"]
        if pmcid and re.fullmatch(r"PMC\d+", pmcid):
            xml_path = RAW / (pmcid + "_fullText.xml")
            provenance = get("https://www.ebi.ac.uk/europepmc/webservices/rest/" + pmcid + "/fullTextXML", xml_path, 5_000_000)
            paper["fulltext_provenance"] = provenance
            if not provenance["error"]:
                article = ET.fromstring(xml_path.read_bytes())
                for element in article.iter("supplementary-material"):
                    hrefs = [v for child in element.iter() for k, v in child.attrib.items() if k.endswith("href")]
                    paper["supplementary_materials"].append({"text": " ".join(element.itertext()).strip()[:1200], "hrefs": hrefs})
                for element in article.iter("sec"):
                    title = element.find("title")
                    if title is not None and re.search(r"data avail|availability of data", " ".join(title.itertext()), re.I):
                        paper["data_availability"].append(" ".join(element.itertext()).strip())
                for element in article.iter("p"):
                    text = " ".join(element.itertext()).strip()
                    if re.search(r"cystatin|creatinine|estimated glomerular|eGFR|GSE\d+|grip strength|leg extension", text, re.I):
                        paper["clinical_excerpts"].append(text[:2500])
        papers.append(paper)
        print(json.dumps({"paper": paper["pmid"], "pmcid": pmcid, "supplements": len(paper["supplementary_materials"]), "fulltext_error": paper.get("fulltext_provenance", {}).get("error")}), flush=True)
    return papers


def main() -> None:
    sys.stdout.reconfigure(encoding="utf-8")
    DATA.mkdir(parents=True, exist_ok=True)
    downloads = []
    for accession, urls in FILES.items():
        for url in urls:
            path = DATA / accession / Path(urlparse(url).path).name
            result = {"accession": accession, **get(url, path)}
            if not result["error"] and ".R.gz" not in path.name:
                result["matrix_preview"] = matrix_preview(path)
            downloads.append(result)
            preview = result.get("matrix_preview", {})
            print(json.dumps({"accession": accession, "file": path.name, "bytes": result.get("bytes"),
                              "data_rows": preview.get("data_rows"), "columns": len(preview.get("header") or []),
                              "header_preview": (preview.get("header") or [])[:8], "error": result.get("error")}), flush=True)
    (OUT / "download_manifest_and_matrix_audit.json").write_text(json.dumps(downloads, ensure_ascii=False, indent=2), encoding="utf-8")
    (OUT / "clinical_metadata_distributions.json").write_text(json.dumps(summarize_samples(), ensure_ascii=False, indent=2), encoding="utf-8")
    if "--matrix-only" not in sys.argv:
        papers = paper_audit()
        (OUT / "candidate_papers_clinical_audit.json").write_text(json.dumps(papers, ensure_ascii=False, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
