"""Donor-level age analysis of the 17-donor human muscle atlas.

The independent unit is donor. Cell and nucleus libraries are summed within
donor x annotation before normalization. This script tests the fixed OXPHOS
genes/module and DoRothEA activities; it does not claim disease or disability
mechanism because the atlas has no functional phenotype.
"""
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

import decoupler as dc
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests

import os
ROOT = Path(os.environ["REPRO_WORK"]) / "analysis"
IN = ROOT / "single_cell_depth_20261008"
DATA = IN / "data"
OUT = IN / "age_tf"
TABLES = OUT / "tables"
AUDIT = OUT / "audit"
REPORTS = OUT / "reports"
GENES = ["COX5A", "MDH2", "SDHA", "GOT2"]
MODULE = GENES
MAJOR = ["MF-I", "MF-II", "MuSC", "FB", "SMC", "CapEC", "VenEC", "Macrophage"]


def fdr(p):
    p = np.asarray(p, dtype=float)
    out = np.full(len(p), np.nan)
    ok = np.isfinite(p)
    if ok.any():
        out[ok] = multipletests(p[ok], method="fdr_bh")[1]
    return out


def fit(y, meta, x="age_numeric"):
    frame = pd.concat([pd.Series(y, index=meta.index, name="y"), meta], axis=1).dropna(subset=["y", x, "sex"])
    if len(frame) < 6 or frame[x].nunique() < 3:
        return {"n_donors": len(frame), "beta_age": np.nan, "ci_low": np.nan, "ci_high": np.nan, "p": np.nan, "rho": np.nan, "rho_p": np.nan, "status": "insufficient"}
    frame["female"] = frame["sex"].astype(str).str.lower().eq("female").astype(float)
    X = sm.add_constant(frame[[x, "female"]].astype(float), has_constant="add")
    model = sm.OLS(frame.y.astype(float), X).fit(cov_type="HC3")
    ci = model.conf_int().loc[x]
    rho, rp = spearmanr(frame[x], frame.y)
    return {"n_donors": len(frame), "beta_age": model.params[x], "ci_low": ci.iloc[0], "ci_high": ci.iloc[1], "p": model.pvalues[x], "rho": rho, "rho_p": rp, "status": "ok"}


def main():
    for p in [TABLES, AUDIT, REPORTS]: p.mkdir(parents=True, exist_ok=True)
    meta = pd.read_csv(DATA / "muscle17_pseudobulk_meta.csv", index_col=0)
    arr = np.load(DATA / "muscle17_pseudobulk_counts.npz", allow_pickle=True)
    genes = pd.Series(arr["gene_symbols"].astype(str))
    counts = arr["counts"]
    # Sum cell/nucleus libraries within donor x annotation. This removes assay
    # as a pseudo-replicate while preserving cell-type-specific profiles.
    meta["age_numeric"] = pd.to_numeric(meta["age_numeric"], errors="coerce")
    meta["sex"] = meta["sex"].astype(str)
    meta["cell_group"] = meta["cell_group"].astype(str)
    keep = meta["n_cells"] >= 30
    meta = meta.loc[keep].copy(); counts = counts[keep.to_numpy()]
    donor_group = meta["donor_id"].astype(str) + "|" + meta["cell_group"].astype(str)
    codes, labels = pd.factorize(donor_group)
    summed = np.zeros((len(labels), counts.shape[1]), dtype=np.float64)
    for i, code in enumerate(codes): summed[code] += counts[i]
    pmeta = meta.groupby(donor_group, observed=True).first().reindex(labels)
    pmeta["donor_group"] = labels
    pmeta.index = labels
    pmeta["library_size"] = summed.sum(axis=1)
    logcpm = np.log2((summed + 0.5) / (summed.sum(axis=1, keepdims=True) + 1.0) * 1e6)
    def col(g):
        hits = np.flatnonzero(genes.str.upper().eq(g.upper()).to_numpy())
        return int(hits[0]) if len(hits) else None
    idx = {g: col(g) for g in GENES}
    missing = [g for g, j in idx.items() if j is None]
    if missing: raise RuntimeError(f"Missing fixed candidate genes: {missing}")
    expression = pd.DataFrame({g: logcpm[:, j] for g, j in idx.items()}, index=pmeta.index)
    expression["OXPHOS_core"] = expression[GENES].apply(lambda x: (x - x.mean()) / x.std(ddof=1) if x.std(ddof=1) > 0 else x * 0).mean(axis=1)
    gene_rows=[]; module_rows=[]
    for group, frame in pmeta.groupby("cell_group", observed=True):
        ii = frame.index
        for g in GENES:
            result=fit(expression.loc[ii,g], frame)
            gene_rows.append({"cell_group":group,"feature":g,**result})
        result=fit(expression.loc[ii,"OXPHOS_core"], frame)
        module_rows.append({"cell_group":group,"feature":"OXPHOS_core",**result})
    gene_df=pd.DataFrame(gene_rows); module_df=pd.DataFrame(module_rows)
    gene_df["fdr_all_celltypes"] = fdr(gene_df["p"]); module_df["fdr_all_celltypes"] = fdr(module_df["p"])
    gene_df.to_csv(TABLES/"candidate_age_associations_by_celltype.csv",index=False)
    module_df.to_csv(TABLES/"oxphos_module_age_associations_by_celltype.csv",index=False)
    # TF activity is computed on donor x cell type pseudobulk only for major
    # classes, with at least six donors. This is a prior-based activity score.
    major = pmeta["cell_group"].isin(MAJOR) & pmeta["age_numeric"].notna()
    matrix = logcpm[major]
    matrix = pd.DataFrame(matrix, index=pmeta.index[major], columns=genes.str.upper())
    # Multiple Ensembl features can map to the same approved symbol; collapse
    # them before passing the matrix to a prior-network scorer.
    matrix = matrix.T.groupby(level=0, sort=True).mean().T
    dor = pd.read_csv(ROOT / "input/dorothea_network.csv")
    dor = dor.rename(columns={"source_genesymbol":"source","target_genesymbol":"target"})
    dor["source"] = dor.source.astype(str).str.upper(); dor["target"] = dor.target.astype(str).str.upper()
    dor = dor[dor.target.isin(matrix.columns)]
    scores, _ = dc.mt.ulm(matrix, dor[["source","target","weight"]], tmin=5, raw=False, verbose=False)
    tf_rows=[]
    for source in scores.columns:
        for group, frame in pmeta.loc[major].groupby("cell_group", observed=True):
            ii=frame.index.intersection(scores.index)
            result=fit(scores.loc[ii,source],frame.loc[ii])
            tf_rows.append({"cell_group":group,"regulator":source,"n_targets":int(dor.loc[dor.source.eq(source),"target"].nunique()),**result})
    tf_df=pd.DataFrame(tf_rows)
    tf_df["fdr_within_regulator_screen"] = tf_df.groupby("cell_group",observed=True)["p"].transform(lambda x: fdr(x))
    tf_df.to_csv(TABLES/"dorothea_tf_age_associations_by_celltype.csv",index=False)
    write = {
        "generated_utc":datetime.now(timezone.utc).isoformat(),"asset":"Human muscle cells and nuclei",
        "dataset_id":"b3c9904b-2ffd-4038-8620-ee1c5b072e9c","donors_total":int(meta.donor_id.nunique()),
        "known_age_donors":int(meta.loc[meta.age_numeric.notna(),"donor_id"].nunique()),
        "independent_unit":"donor; cell and nucleus assays summed within donor x annotation",
        "minimum_cells_per_assay_group":30,"features":GENES+['OXPHOS_core'],
        "cell_groups_tested":sorted(module_df.cell_group.unique().tolist()),
        "major_tf_cell_groups":MAJOR,
        "tf_method":"DoRothEA A/B/C + decoupler ULM on donor x cell-group log2 CPM",
        "age_model":"activity or expression ~ numeric age + sex, HC3 covariance",
        "interpretation_limit":"Normal intercostal muscle atlas; age trends do not establish disability, renal discordance, disease regulation, or mediation.",
        "assay_scope":"cells and nuclei were summed; remaining batch and tissue processing confounding requires sensitivity analysis."
    }
    (AUDIT/"age_tf_analysis_specification.json").write_text(json.dumps(write,ensure_ascii=False,indent=2),encoding="utf-8")
    report=["# 多供者肌肉单细胞：供者级年龄与 TF 活性分析","","## 统计单位和模型","","将细胞与细胞核测序按 donor × annotation 汇总，供者作为独立样本；每个特征用 `表达或活性 ~ 年龄 + 性别` 的 HC3 回归。只把四个预先固定的 OXPHOS 候选基因和核心模块作为主要屏幕，DoRothEA/ULM 作为上游调控探索。","","## OXPHOS 模块年龄关联","","| 细胞群 | donor 数 | 年龄 beta | 95% CI | P | FDR |","|---|---:|---:|---:|---:|---:|"]
    for r in module_df.sort_values("p",na_position="last").itertuples(): report.append(f"| {r.cell_group} | {r.n_donors} | {r.beta_age:.3f} | [{r.ci_low:.3f}, {r.ci_high:.3f}] | {r.p:.3g} | {r.fdr_all_celltypes:.3g} |" if np.isfinite(r.p) else f"| {r.cell_group} | {r.n_donors} | NA | NA | NA | NA |")
    report += ["","## 解释边界","","该资源包含多个年龄阶段和 16 名可数值化年龄供者，能够把原单供者定位升级为供者级年龄趋势证据；但没有 ADL、衰弱、肾功能或疾病对比，不能称为当前临床主线的直接机制验证。拟时序若用于 MuSC-肌纤维分化，只能是状态连续性探索；年龄回归和供者级伪 bulk 是更强的独立证据层。","","输出文件：`candidate_age_associations_by_celltype.csv`、`oxphos_module_age_associations_by_celltype.csv`、`dorothea_tf_age_associations_by_celltype.csv`。"]
    (REPORTS/"multidonor_age_tf_analysis_zh.md").write_text("\n".join(report)+"\n",encoding="utf-8")
    print("MULTIDONOR_AGE_TF_COMPLETE")
    print(module_df.sort_values("p",na_position="last").head(20).to_string(index=False))

if __name__ == "__main__": main()
