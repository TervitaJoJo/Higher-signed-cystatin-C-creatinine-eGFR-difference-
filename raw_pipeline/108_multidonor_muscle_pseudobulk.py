"""Audit frozen public H5AD assets and aggregate raw counts by donor and cell type.

Remote metadata use AnnData's HDF5 reader; expression is processed locally in
CSR blocks to avoid loading both atlas matrices into memory simultaneously.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import time
from datetime import datetime, timezone
from pathlib import Path

import anndata as ad
import fsspec
import h5py
import numpy as np
import pandas as pd
import requests
from scipy import sparse

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "single_cell_depth_20261008"
DATA = OUT / "data"
TABLES = OUT / "tables"
AUDIT = OUT / "audit"
ASSETS = {
    "muscle17": ("b3c9904b-2ffd-4038-8620-ee1c5b072e9c",
                 "https://datasets.cellxgene.cziscience.com/5fd93c55-8047-48fb-905c-6ef2f7767483.h5ad"),
    "muscle32": ("15d374d6-0dfd-4d8e-ade7-81f73dc921ee",
                 "https://datasets.cellxgene.cziscience.com/1cb67ec8-c9e3-4d4a-ba23-a1a8ef3f8450.h5ad"),
}


def age_number(row: pd.Series) -> float:
    for key in ["age", "Age", "age_years", "Age_bin", "development_stage"]:
        text = str(row.get(key, ""))
        m = re.search(r"(\d+(?:\.\d+)?)\s*[- ]?(?:year|years|y/o)", text.lower())
        if m:
            return float(m.group(1))
        m = re.fullmatch(r"\s*(\d+(?:\.\d+)?)\s*", text)
        if m:
            return float(m.group(1))
        m = re.search(r"(\d+)\s*[-_to]+\s*(\d+)", text)
        if m:
            return (float(m.group(1)) + float(m.group(2))) / 2
        decades = {"third": 25, "fourth": 35, "fifth": 45, "sixth": 55,
                   "seventh": 65, "eighth": 75, "ninth": 85}
        for word, middle in decades.items():
            if word + " decade" in text.lower():
                return float(middle)
    return np.nan


def write_json(path: Path, obj: dict) -> None:
    path.write_text(json.dumps(obj, indent=2, default=str), encoding="utf-8")


def inspect(handle: h5py.File, name: str) -> tuple[pd.DataFrame, pd.DataFrame, str]:
    obs = ad.io.read_elem(handle["obs"])
    var = ad.io.read_elem(handle["var"])
    if "raw" in handle and "X" in handle["raw"]:
        raw = "raw/X"
        rawvar = ad.io.read_elem(handle["raw/var"])
    elif "layers" in handle and "counts" in handle["layers"]:
        raw = "layers/counts"
        rawvar = var
    else:
        raw = "X"
        rawvar = var
    print(name, "cells", len(obs), "donors", obs.donor_id.nunique(), "matrix", raw, flush=True)
    obs.to_csv(DATA / f"{name}_obs.csv.gz", compression="gzip")
    rawvar.to_csv(DATA / f"{name}_var.csv.gz", compression="gzip")
    obs["age_numeric"] = obs.apply(age_number, axis=1)
    donor_columns = [c for c in ["donor_id", "age_numeric", "Age_bin", "Age_group", "age", "Age", "sex",
        "assay", "suspension_type", "tissue", "batch", "SampleID", "sample", "development_stage"] if c in obs]
    donors = obs.groupby(donor_columns, observed=True, dropna=False).size().rename("n_cells").reset_index()
    donors.to_csv(TABLES / f"{name}_donor_metadata.csv", index=False)
    labels = [c for c in ["annotation_level0", "annotation_level1", "annotation_level2", "cell_type", "cell_type_original"] if c in obs]
    for col in labels:
        obs.groupby(["donor_id", col, "assay"], observed=True, dropna=False).size().rename("n_cells").reset_index().to_csv(
            TABLES / f"{name}_{col}_counts.csv", index=False)
    obj = handle[raw]
    vals = obj["data"][:10000] if isinstance(obj, h5py.Group) else obj[:100, :100].ravel()
    encoding = dict(obj.attrs)
    audit = {"asset": ASSETS[name], "n_cells": len(obs), "donors": obs.donor_id.nunique(),
        "obs_columns": list(obs.columns), "matrix": raw, "matrix_attrs": encoding,
        "var_columns": list(rawvar.columns), "first_nonzero_integer_fraction": float(np.mean(np.isclose(vals, np.rint(vals)))),
        "age_known_donors": obs.loc[obs.age_numeric.notna(), "donor_id"].nunique(),
        "age_range": [float(obs.age_numeric.min()), float(obs.age_numeric.max())],
        "annotations": {c: obs[c].value_counts().to_dict() for c in labels},
        "generated_utc": datetime.now(timezone.utc).isoformat()}
    write_json(AUDIT / f"{name}_metadata_audit.json", audit)
    print(donors.to_string(index=False), flush=True)
    # Keep console output ASCII-safe on Windows; full labels remain in CSV/JSON.
    print({c: len(obs[c].value_counts()) for c in labels}, flush=True)
    return obs, rawvar, raw


def download(name: str) -> Path:
    url = ASSETS[name][1]
    path = DATA / f"{name}.h5ad"
    if path.exists():
        return path
    partial = path.with_suffix(".partial")
    session = requests.Session()
    for attempt in range(5):
        done = partial.stat().st_size if partial.exists() else 0
        headers = {"Range": f"bytes={done}-"} if done else {}
        try:
            response = session.get(url, headers=headers, stream=True, timeout=(30, 90))
            response.raise_for_status()
            append = done > 0 and response.status_code == 206
            if done and not append:
                done = 0
            total = done + int(response.headers.get("Content-Length", 0))
            previous = time.monotonic()
            with partial.open("ab" if append else "wb") as dst:
                for block in response.iter_content(chunk_size=4 * 1024**2):
                    dst.write(block)
                    done += len(block)
                    if time.monotonic() - previous > 15:
                        print(f"{name} download {done / 1e6:.0f}/{total / 1e6:.0f} MB", flush=True)
                        previous = time.monotonic()
            if total and done != total:
                raise IOError("Incomplete asset transfer")
            with h5py.File(partial, "r") as check:
                assert "obs" in check and "var" in check
            partial.replace(path)
            sha = hashlib.file_digest(path.open("rb"), "sha256").hexdigest()
            write_json(AUDIT / f"{name}_download.json", {"url": url, "bytes": done, "sha256": sha,
                "downloaded_utc": datetime.now(timezone.utc).isoformat()})
            return path
        except (requests.RequestException, IOError) as exc:
            print(name, "retry", attempt + 1, str(exc), flush=True)
    raise RuntimeError(f"Download failed: {name}")


def aggregate(name: str, path: Path) -> None:
    with h5py.File(path, "r") as handle:
        obs, var, raw = inspect(handle, name)
        obj = handle[raw]
        assert obj.attrs["encoding-type"] == "csr_matrix", "Expected CSR raw-count matrix"
        if name == "muscle17":
            cellkey = "annotation_level1"
        else:
            cellkey = "cell_type"
        keys = ["donor_id", cellkey, "assay"]
        columns = keys + [c for c in ["age_numeric", "Age_group", "Age_bin", "sex", "tissue", "batch", "SampleID", "suspension_type"] if c in obs and c not in keys]
        valid = obs["is_primary_data"].astype(bool).to_numpy() if "is_primary_data" in obs else np.ones(len(obs), dtype=bool)
        group_index = pd.MultiIndex.from_frame(obs[keys].astype(str))
        codes, groups = pd.factorize(group_index)
        ng, nv = len(groups), len(var)
        counts = np.zeros((ng, nv), dtype=np.float64)
        detection = np.zeros((ng, nv), dtype=np.int32)
        ncells = np.zeros(ng, dtype=np.int64)
        library = np.zeros(len(obs))
        ptr = obj["indptr"][:]
        integer_bad = 0
        for start in range(0, len(obs), 4000):
            stop = min(len(obs), start + 4000)
            left, right = int(ptr[start]), int(ptr[stop])
            data = obj["data"][left:right]
            integer_bad += int(np.sum(~np.isclose(data, np.rint(data))))
            assert np.all(data >= 0) and np.all(np.isfinite(data))
            block = sparse.csr_matrix((data, obj["indices"][left:right], ptr[start:stop+1] - left), shape=(stop-start, nv))
            # Aggregate with sparse incidence multiplication; each donor remains one replicate.
            use = valid[start:stop]
            incidence = sparse.csr_matrix((use.astype(float), (codes[start:stop], np.arange(stop-start))), shape=(ng, stop-start))
            counts += (incidence @ block).toarray()
            detected = block.copy()
            detected.data = (detected.data > 0).astype(np.int32)
            detection += (incidence @ detected).toarray().astype(np.int32)
            ncells += np.bincount(codes[start:stop][use], minlength=ng)
            library[start:stop] = np.asarray(block.sum(axis=1)).ravel()
            if start % 20000 == 0:
                print(name, "aggregated", stop, "/", len(obs), flush=True)
        assert integer_bad == 0, f"{integer_bad} non-integer values: matrix is not raw counts"
        meta = obs[columns].assign(pb_index=codes).groupby("pb_index", observed=True).first().reindex(range(ng))
        meta = meta.rename(columns={cellkey: "cell_group"})
        meta["n_cells"] = ncells
        meta["library_size"] = counts.sum(axis=1)
        meta.index.name = "pb_index"
        meta.to_csv(DATA / f"{name}_pseudobulk_meta.csv")
        np.savez_compressed(DATA / f"{name}_pseudobulk_counts.npz", counts=counts, detection=detection,
                            gene_symbols=var["feature_name"].astype(str).to_numpy(), feature_ids=var.index.astype(str).to_numpy())
        candidate_columns = [int(np.flatnonzero(var.feature_name.astype(str).to_numpy() == g)[0]) for g in ["COX5A", "MDH2", "SDHA", "GOT2"]]
        qc = obs[["donor_id", cellkey, "assay"]].copy()
        qc["library_size"] = library
        qc.to_csv(DATA / f"{name}_cell_library_sizes.csv.gz", compression="gzip")
        write_json(AUDIT / f"{name}_pseudobulk_audit.json", {"groups": ng, "raw_count_path": raw,
            "primary_cells": int(valid.sum()), "gene_count": nv, "noninteger_values": integer_bad,
            "groups_by": keys, "candidate_columns": candidate_columns, "count_sum": float(counts.sum()),
            "minimum_cells_for_inference": 30, "no_new_cell_filter": "Atlas author QC retained; raw counts verified"})


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=["metadata", "aggregate"], required=True)
    parser.add_argument("--datasets", nargs="+", default=list(ASSETS))
    args = parser.parse_args()
    for folder in [DATA, TABLES, AUDIT]:
        folder.mkdir(parents=True, exist_ok=True)
    for name in args.datasets:
        if args.mode == "metadata":
            local = DATA / f"{name}.h5ad"
            if local.exists():
                with h5py.File(local, "r") as handle:
                    inspect(handle, name)
            else:
                with fsspec.open(ASSETS[name][1], block_size=4*1024**2).open() as remote:
                    with h5py.File(remote, "r") as handle:
                        inspect(handle, name)
        else:
            aggregate(name, download(name))


if __name__ == "__main__":
    main()
