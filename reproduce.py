'''Portable reproduction from processed observations or locked aggregate sources.'''
from __future__ import annotations
import argparse
import hashlib
import io
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time

sys.dont_write_bytecode = True
P = Path(__file__).resolve().parent
RESULTS = P / "results"

def read_tables(name):
    return json.loads((P/"data"/name).read_text(encoding="utf-8"))

def verify_source_checksums():
    manifest=P/'checksums.sha256'
    if not manifest.exists():return
    for line in manifest.read_text(encoding='utf-8').splitlines():
        digest,name=line.split('  ',1)
        if not name.startswith(('data/','scripts/','raw_pipeline/')):continue
        file=P/name
        if not file.is_file() or hashlib.sha256(file.read_bytes()).hexdigest()!=digest:
            raise RuntimeError('Source integrity check failed: '+name)

def write_sources(work):
    a=work/"analysis"
    for name,text in read_tables("reference_tables.json").items():
        target=a/name;target.parent.mkdir(parents=True,exist_ok=True)
        target.write_text(text,encoding="utf-8")
    source=a/"audit/bmc_figure_restructure_20261011/source_data"
    source.mkdir(parents=True,exist_ok=True)
    for name,text in read_tables("figure_source_data.json").items():
        (source/name).write_text(text,encoding="utf-8")
    for folder in ["data_internal","tables","audit","logs","reports","input","code"]:
        (a/folder).mkdir(exist_ok=True)
    for target in [work/"BMC geriatrics/V1/figures/panels",work/"BMC geriatrics/V1/supplementary_figures"]:
        target.mkdir(parents=True,exist_ok=True)
    for f in (P/"scripts").glob("*.R"): shutil.copyfile(f,a/"code"/f.name)
    shutil.copyfile(P/"data/bulk_expression_inputs.rds",a/"input/bulk_expression_inputs.rds")
    resource=a/"transcriptome_analysis/resources";resource.mkdir(parents=True,exist_ok=True)
    shutil.copyfile(P/"data/MSigDB_Hallmark_2020.gmt",resource/"MSigDB_Hallmark_2020.gmt")
    data=a/"single_cell_depth_20261008/data";data.mkdir(parents=True,exist_ok=True)
    shutil.copyfile(P/"data/muscle_pseudobulk_counts.npz",data/"muscle17_pseudobulk_counts.npz")
    shutil.copyfile(P/"data/muscle_pseudobulk_metadata.csv",data/"muscle17_pseudobulk_meta.csv")
    target=a/"input"
    shutil.copyfile(P/"data/dorothea_network.csv",target/"dorothea_network.csv")

def rcall(work,script,*args):
    rscript=os.environ.get("RSCRIPT",shutil.which("Rscript") or "")
    if not rscript: raise RuntimeError("Rscript was not found. Install R or set RSCRIPT.")
    script_path=work/"analysis/code"/script
    expr='suppressWarnings(Sys.setlocale("LC_CTYPE","English_United States.utf8")); .review_args<-commandArgs(TRUE); commandArgs<-function(trailingOnly=FALSE) if(trailingOnly) .review_args[-1] else base::commandArgs(FALSE); if(nzchar(Sys.getenv("REVIEW_R_LIB"))) .libPaths(c(Sys.getenv("REVIEW_R_LIB"),.libPaths())); source(.review_args[1],chdir=FALSE)'
    return command([rscript,"-e",expr,str(script_path),str(work/"analysis"),*map(str,args)],work,script)

def command(args,work,label):
    start=time.monotonic();print("Running",label,flush=True)
    log=work/(label.replace(".R","").replace(".py","")+".log")
    with log.open("w",encoding="utf-8") as f:
        process=subprocess.run(args,cwd=work,env=os.environ,stdout=f,stderr=subprocess.STDOUT)
    if process.returncode:
        print("\n".join(log.read_text(encoding="utf-8",errors="replace").splitlines()[-18:]),file=sys.stderr)
        raise RuntimeError(label+" failed; inspect the temporary working directory.")
    print(f"Completed {label} ({time.monotonic()-start:.1f}s)",flush=True)

def compare_tables(work,mode):
    import numpy as np
    import pandas as pd
    rows=[]
    for name,text in read_tables("reference_tables.json").items():
        is_clinical=not name.startswith(("muscle_","single_cell_"))
        expected=mode=="full" or (mode=="public" and not is_clinical)
        if not expected: continue
        path=work/"analysis"/name
        if not path.exists():
            rows.append(dict(table=name,status="missing",numeric_cells=0,max_abs_error=None));continue
        old=pd.read_csv(io.StringIO(text));new=pd.read_csv(path)
        if old.shape!=new.shape or list(old.columns)!=list(new.columns):
            rows.append(dict(table=name,status="shape_or_columns_differ",numeric_cells=0,max_abs_error=None));continue
        error=0.;cells=0;passed=True
        for col in old:
            if pd.api.types.is_numeric_dtype(old[col]):
                a=old[col].to_numpy(float);b=pd.to_numeric(new[col],errors="coerce").to_numpy(float)
                if not np.allclose(a,b,rtol=1e-7,atol=1e-9,equal_nan=True): passed=False
                ok=np.isfinite(a)&np.isfinite(b)
                if ok.any():error=max(error,float(np.max(np.abs(a[ok]-b[ok]))));cells+=int(ok.sum())
            elif not old[col].fillna("").astype(str).equals(new[col].fillna("").astype(str)):
                passed=False
        rows.append(dict(table=name,status="passed" if passed else "different",numeric_cells=cells,max_abs_error=error))
    return rows

def draw_all(work):
    import matplotlib
    matplotlib.use("Agg")
    sys.path.insert(0,str(P/"scripts"))
    import plot_base as d
    import plot_clinical as c
    import plot_matrix as m
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    c.COHORTS=d.COHORTS
    original=d.save
    def clinical_save(fig,name):
        for ax in fig.axes:
            for text in list(ax.texts):
                if text.get_text().startswith('Public skeletal muscle transcriptomes:'): text.remove()
        original(fig,name)
    d.save=clinical_save;d.panel1a();d.save=original
    d.panel1b();d.panel1c()
    d.assemble('Figure1_design_risk_set_and_states',180,178,
        [('Figure1A_risk_set_flow',0,0,180,95),('Figure1B_observed_states',0,95,90,83),('Figure1C_exposure_distribution',90,95,90,83)])
    c.panel2a();d.panel2b();d.panel2c()
    d.assemble('Figure2_recovery_association_shape_and_probabilities',180,160,
        [('Figure2A_clinical_associations',0,0,180,72),('Figure2B_existing_recovery_splines',0,72,90,88),('Figure2C_probability_contrasts',90,72,90,88)])
    c.panel3a()
    rows=[c.model(cohort,analysis,sample) for cohort in ['HRS2016VBS','CHARLS2011']
        for analysis,sample in [('combined_base','bun_hgb_alb' if cohort=='HRS2016VBS' else 'bun_hgb'),
        ('combined_bun_hgb_alb' if cohort=='HRS2016VBS' else 'combined_bun_hgb','bun_hgb_alb' if cohort=='HRS2016VBS' else 'bun_hgb')]]
    m.draw_panel(rows)
    d.assemble('Figure3_clinical_sensitivity_comparisons',180,178,
        [('Figure3A_filtration_adjustment_paths',0,0,180,88),('Figure3B_matched_case_estimate_matrix',0,88,180,90)])
    d.supplementary_shape();d.supplementary_nonlinearity()
    d.panel3a();d.panel3b();d.panel3c()
    d.assemble('FigureS6_external_muscle_OXPHOS_context',180,154,
        [('Figure3A_OXPHOS_family_forest',0,0,180,68),('Figure3B_candidate_directions',0,68,90,86),('Figure3C_leave_one_family_out',90,68,90,86)],d.V1/'supplementary_figures')
    def forest(frame,stem,labelcol,estimate,lower,upper,title,xlabel):
        z=frame.reset_index(drop=True);fig,ax=plt.subplots(figsize=(7.1,max(2.5,len(z)*.36+1)))
        y=np.arange(len(z));colors=[d.HRS if str(x).startswith('HRS') else d.CHARLS for x in z.get('cohort',pd.Series(['']*len(z)))]
        for j,row in z.iterrows():ax.errorbar(row[estimate],y[j],xerr=[[row[estimate]-row[lower]],[row[upper]-row[estimate]]],fmt='o',color=colors[j],capsize=2)
        ax.axvline(1 if estimate=='odds_ratio' else 0,color='#999999',ls='--',lw=.8)
        ax.set_yticks(y,z[labelcol]);ax.invert_yaxis();ax.set_xlabel(xlabel);ax.set_title(title,fontsize=10);d.clean(ax);fig.tight_layout()
        for ext in ['pdf','png']:fig.savefig(d.V1/'supplementary_figures'/f'{stem}.{ext}',dpi=300)
        plt.close(fig)
    age=pd.read_csv(work/'analysis/mortality_heterogeneity_20261008/tables/mortality_stratified_results.csv')
    age=age[(age.cohort=='HRS2016VBS') & age.subgroup.str.startswith('age_group')].copy()
    overall=pd.read_csv(work/'analysis/tables/final_competing_transition_results.csv')
    overall=overall[(overall.cohort=='CHARLS2011') & (overall.outcome=='Death')].rename(columns={'ci_lower':'ci_low','ci_upper':'ci_high'})
    overall['subgroup']='All ages'
    age=pd.concat([age,overall],ignore_index=True)
    age['display']=age.cohort.astype(str)+' '+age.subgroup.astype(str).str.replace('age_group=interval_start_age_','',regex=False)
    forest(age,'FigureS2_mortality_age_heterogeneity','display','odds_ratio','ci_low','ci_high','Exploratory mortality associations','Adjacent-wave mortality OR')
    module=pd.read_csv(work/'analysis/single_cell_depth_20261008/age_tf/tables/oxphos_module_age_associations_by_celltype.csv')
    module=module.set_index('cell_group').loc[['MF-I','MF-II','MF-Isc(fg)','MuSC','FB','MTJ']].reset_index()
    forest(module,'FigureS4_single_cell_OXPHOS','cell_group','beta_age','ci_low','ci_high','Donor-level muscle OXPHOS age associations','Age coefficient per year')
    lr=pd.read_csv(work/'analysis/single_cell_depth_20261008/communication/tables/canonical_ligand_receptor_age_screen.csv')
    lr=lr.sort_values('fdr_screen').drop_duplicates(['ligand','receptor']).copy()
    lr['pair']=lr.ligand+' - '+lr.receptor
    forest(lr,'FigureS5_ligand_receptor_screen','pair','beta_age','ci_low','ci_high','Exploratory expression co-availability','Age coefficient per year')
    rendered=[]
    for folder in ['figures','supplementary_figures']:
        for f in sorted((d.V1/folder).glob('*')):
            if f.is_file() and f.suffix in ['.pdf','.png']:
                target=RESULTS/f.name;shutil.copyfile(f,target);rendered.append(target.name)
    return rendered

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode',choices=['figures','full','public'],default='figures')
    parser.add_argument('--work-dir',type=Path)
    parser.add_argument('--keep-work',action='store_true')
    args=parser.parse_args()
    verify_source_checksums()
    RESULTS.mkdir(exist_ok=True)
    work=args.work_dir.resolve() if args.work_dir else Path(tempfile.mkdtemp(prefix='adl-review-'))
    if work==P or P in work.parents: raise ValueError('Work directory must be outside the package to enforce the 100-file limit.')
    work.mkdir(parents=True,exist_ok=True)
    os.environ['REPRO_WORK']=str(work);os.environ['PYTHONDONTWRITEBYTECODE']='1'
    default_r_library=Path.home()/'R/review-reproduction-library'
    if 'REVIEW_R_LIB' not in os.environ and default_r_library.is_dir():
        os.environ['REVIEW_R_LIB']=str(default_r_library)
    success=False
    try:
        write_sources(work)
        a=work/'analysis'
        for name in read_tables('reference_tables.json'):
            if args.mode=='full' or (args.mode=='public' and name.startswith(('muscle_','single_cell_'))):
                (a/name).unlink()
        if args.mode=='full':
            for name in ['clinical_panels.rds','renal_intervals.rds','selection_intervals.rds']:
                src=P/'local_data'/name
                if not src.is_file(): raise FileNotFoundError('Full model reproduction requires lawfully obtained processed clinical inputs; see DATA_ACCESS.md.')
                shutil.copyfile(src,a/'input'/name)
            for step in ['prepare_clinical.R','102_clinical_nonlinear_interaction_competing_risk.R','materialise_primary.R','73_complete_case_selection.R','106_mortality_heterogeneity.R',
                '113_hrs_primary_evalue_and_framing.R','119_renal_strata_discordance.R','124_renal_extension.R',
                '130_primary_baseline_characteristics.R','105_incremental_prediction_gain.R','110_prediction_gain_design_validation.R',
                'export_clinical_visuals.R']:
                rcall(work,step)
        if args.mode in ['full','public']:
            for step in ['98_muscle_oxphos_reaudit.R','101_muscle_oxphos_gene_stability.R','104_project_level_transcriptome_validation.R']:
                rcall(work,step)
            for step in ['single_cell_oxphos.py','single_cell_lr.py']:
                command([sys.executable,str(P/'scripts'/step)],work,step)
        verification=compare_tables(work,args.mode)
        (RESULTS/'verification.json').write_text(json.dumps({'mode':args.mode,'comparison_tolerance':{'rtol':1e-7,'atol':1e-9},'tables':verification},indent=2),encoding='utf-8')
        failures=[r for r in verification if r['status']!='passed']
        if failures: raise RuntimeError(f'{len(failures)} result tables differ; see results/verification.json. Locked references are unchanged.')
        figures=draw_all(work)
        reproduced={}
        for name in read_tables('reference_tables.json'):
            path=a/name
            if path.exists(): reproduced[name]=path.read_text(encoding='utf-8-sig')
        (RESULTS/'reproduced_tables.json').write_text(json.dumps(reproduced,indent=1),encoding='utf-8')
        summary={'status':'passed','mode':args.mode,'tables_compared':len(verification),
            'numeric_cells_compared':sum(r['numeric_cells'] for r in verification),
            'figure_files':figures,'model_refitting':args.mode in ['full','public'],
            'clinical_models_refitted':args.mode=='full','original_manuscript_modified':False}
        (RESULTS/'run_summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
        success=True
        count=sum(f.is_file() for f in P.rglob('*'))
        if count>100:raise RuntimeError(f'Package has {count} files; limit is 100.')
        print(json.dumps(summary,indent=2));print('Package files:',count)
    finally:
        if success and not args.keep_work:shutil.rmtree(work)
        else:print('Working directory retained:',work)

if __name__=='__main__':main()
