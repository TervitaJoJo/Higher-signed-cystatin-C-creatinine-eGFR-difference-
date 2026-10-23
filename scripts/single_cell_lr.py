"""Bounded donor-level ligand-receptor screen in the normal muscle atlas.

This is a prespecified hypothesis screen, not a CellChat causal network: each
donor x cell-group logCPM profile is used to score expression of a curated
ligand/receptor pair, and scores are tested against age with donor as unit.
"""
from __future__ import annotations
import json
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests

import os
ROOT=Path(os.environ["REPRO_WORK"])/"analysis"; IN=ROOT/"single_cell_depth_20261008"; DATA=IN/"data"; OUT=IN/"communication"; TABLES=OUT/"tables"; AUDIT=OUT/"audit"; REPORTS=OUT/"reports"
PAIRS=[
    ("IGF1","IGF1R","growth/regeneration"),("VEGFA","KDR","vascular support"),
    ("CXCL12","CXCR4","progenitor niche"),("TGFB1","TGFBR1","fibrosis/remodelling"),
    ("TGFB1","TGFBR2","fibrosis/remodelling"),("FGF2","FGFR1","regeneration"),
    ("HGF","MET","regeneration"),("CCL2","CCR2","myeloid recruitment"),
    ("TNF","TNFRSF1A","inflammation"),("IL6","IL6R","inflammation")]
SENDERS=["FB","SMC","Macrophage","CapEC","VenEC","MF-I","MF-II","MuSC"]
RECEIVERS=["MuSC","MF-I","MF-II","FB","CapEC","VenEC"]

def fdr(p):
    p=np.asarray(p,float); out=np.full(len(p),np.nan); ok=np.isfinite(p)
    if ok.any(): out[ok]=multipletests(p[ok],method="fdr_bh")[1]
    return out

def model(score, meta):
    z=pd.concat([pd.Series(score,index=meta.index,name="score"),meta],axis=1).dropna(subset=["score","age_numeric","sex"])
    if len(z)<6 or z.age_numeric.nunique()<3: return {"n_donors":len(z),"beta_age":np.nan,"ci_low":np.nan,"ci_high":np.nan,"p":np.nan,"rho":np.nan,"rho_p":np.nan}
    z["female"]=z.sex.astype(str).str.lower().eq("female").astype(float)
    X=sm.add_constant(z[["age_numeric","female"]],has_constant="add"); fit=sm.OLS(z.score,X).fit(cov_type="HC3"); ci=fit.conf_int().loc["age_numeric"]; rho,rp=spearmanr(z.age_numeric,z.score)
    return {"n_donors":len(z),"beta_age":fit.params.age_numeric,"ci_low":ci.iloc[0],"ci_high":ci.iloc[1],"p":fit.pvalues.age_numeric,"rho":rho,"rho_p":rp}

def main():
    for p in [TABLES,AUDIT,REPORTS]: p.mkdir(parents=True,exist_ok=True)
    meta=pd.read_csv(DATA/"muscle17_pseudobulk_meta.csv",index_col=0)
    arr=np.load(DATA/"muscle17_pseudobulk_counts.npz",allow_pickle=True); genes=pd.Series(arr["gene_symbols"].astype(str)); counts=arr["counts"]
    meta["age_numeric"]=pd.to_numeric(meta.age_numeric,errors="coerce"); meta["cell_group"]=meta.cell_group.astype(str)
    keep=meta.n_cells>=30; meta=meta.loc[keep].copy(); counts=counts[keep.to_numpy()]
    donor_group=meta.donor_id.astype(str)+"|"+meta.cell_group.astype(str); codes,labels=pd.factorize(donor_group)
    summed=np.zeros((len(labels),counts.shape[1]))
    for i,c in enumerate(codes): summed[c]+=counts[i]
    pmeta=meta.groupby(donor_group,observed=True).first().reindex(labels); pmeta.index=labels; pmeta["library_size"]=summed.sum(axis=1)
    logcpm=np.log2((summed+0.5)/(summed.sum(axis=1,keepdims=True)+1)*1e6); symbols=genes.astype(str).str.upper()
    matrix=pd.DataFrame(logcpm,index=labels,columns=symbols).T.groupby(level=0,sort=True).mean().T
    result=[]
    for ligand,receptor,axis in PAIRS:
      if ligand not in matrix or receptor not in matrix: continue
      for sender in SENDERS:
       for receiver in RECEIVERS:
        smeta=pmeta[pmeta.cell_group.eq(sender)&pmeta.age_numeric.notna()]; rmeta=pmeta[pmeta.cell_group.eq(receiver)&pmeta.age_numeric.notna()]
        common=smeta.index.intersection(rmeta.index)
        if len(common)<6: continue
        # Product of positive logCPM after shifting each gene to zero; this is
        # an interpretable expression co-availability score, not communication.
        l=matrix.loc[common,ligand].to_numpy(); r=matrix.loc[common,receptor].to_numpy(); score=np.sqrt(np.maximum(l-l.min(),0.0)*np.maximum(r-r.min(),0.0))
        result.append({"ligand":ligand,"receptor":receptor,"axis":axis,"sender":sender,"receiver":receiver,**model(score,smeta.loc[common])})
    df=pd.DataFrame(result)
    if not df.empty: df["fdr_screen"]=fdr(df.p)
    df.to_csv(TABLES/"canonical_ligand_receptor_age_screen.csv",index=False)
    audit={"generated_utc":datetime.now(timezone.utc).isoformat(),"independent_unit":"donor","pairs":PAIRS,"senders":SENDERS,"receivers":RECEIVERS,"score":"sqrt(shifted ligand logCPM x shifted receptor logCPM)","model":"score ~ numeric age + sex, HC3 covariance","atlas":"normal intercostal muscle, 17 donors (16 numeric ages)","interpretation_limit":"Expression co-availability is not a measured ligand-receptor interaction, cell communication flux, or causal mechanism; no disease/ADL phenotype is available."}
    (AUDIT/"canonical_lr_age_screen_specification.json").write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding="utf-8")
    report=["# 预设配体-受体年龄屏幕","","该分析只把公开正常肌肉图谱中的配体和受体共表达作为假设筛选。供者是独立单位，细胞/细胞核先汇总到 donor × 细胞群；没有把表达共现写成真实的细胞通讯。","","| 配体 | 受体 | 发送者 | 接收者 | donor 数 | 年龄 beta | 95% CI | P | FDR |","|---|---|---|---|---:|---:|---:|---:|---:|"]
    if not df.empty:
      for r in df.sort_values("p",na_position="last").head(30).itertuples(): report.append(f"| {r.ligand} | {r.receptor} | {r.sender} | {r.receiver} | {r.n_donors} | {r.beta_age:.3f} | [{r.ci_low:.3f}, {r.ci_high:.3f}] | {r.p:.3g} | {r.fdr_screen:.3g} |" if np.isfinite(r.p) else f"| {r.ligand} | {r.receptor} | {r.sender} | {r.receiver} | {r.n_donors} | NA | NA | NA | NA |")
    report += ["","解释：该屏幕只能提供正常肌肉中年龄相关的配体/受体表达共可用性，不能连接到 eGFR discordance 或 ADL 恢复，也不能替代 CellChat 的疾病对照或功能实验。"]
    (REPORTS/"canonical_lr_age_screen_zh.md").write_text("\n".join(report)+"\n",encoding="utf-8")
    print("CANONICAL_LR_AGE_SCREEN_COMPLETE"); print(df.sort_values("p",na_position="last").head(20).to_string(index=False))
if __name__=="__main__": main()
