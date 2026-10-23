from pathlib import Path
import sys
import json
import hashlib
import shutil
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle, FancyArrowPatch, Patch
from matplotlib.ticker import FixedLocator, FixedFormatter
import fitz
from PIL import Image

from nature_style import apply_style

import os
ROOT = Path(os.environ["REPRO_WORK"])
A = ROOT / 'analysis'
V1 = ROOT / 'BMC geriatrics' / 'V1'
AUDIT = ROOT / "analysis/audit/bmc_figure_restructure_20261011"
SOURCE = AUDIT / 'source_data'
PANELS = V1 / 'figures' / 'panels'
PANELS.mkdir(parents=True, exist_ok=True)
apply_style(body=9, tick=8.5, label=8.5)
plt.rcParams.update({'pdf.fonttype': 42, 'ps.fonttype': 42, 'svg.fonttype': 'none',
                     'axes.unicode_minus': True, 'savefig.bbox': None})
HRS = '#D55E00'
CHARLS = '#0072B2'
RECOVERY = '#299D74'
PERSISTENT = '#E3BE47'
DEATH = '#8B9095'
INK = '#242628'
LOWER = '#277A79'
HIGHER = '#C25C83'
COHORTS = [('HRS2016VBS', 'HRS', HRS), ('CHARLS2011', 'CHARLS', CHARLS)]
MM = 1 / 25.4
manifest = {}


def source(path):
    path = Path(path)
    manifest[str(path.relative_to(ROOT))] = hashlib.sha256(path.read_bytes()).hexdigest()
    return pd.read_csv(path)


clinical = source(A / 'tables/final_competing_transition_results.csv')
enh = source(A / 'tables/clinical_enhancement_model_results.csv')
nonlinear = source(A / 'tables/clinical_nonlinear_spline_results.csv')
interactions = source(A / 'tables/clinical_interaction_results.csv')
renal = source(A / 'renal_extension_20261011/tables/renal_extension_model_results.csv')
flow = source(SOURCE / 'risk_set_flow.csv')
histogram = source(SOURCE / 'risk_set_exposure_histogram.csv')
quantiles = source(SOURCE / 'risk_set_exposure_quantiles.csv')
curves = source(SOURCE / 'existing_recovery_spline_curves.csv')
muscle = A / 'muscle_oxphos_reaudit_20261007/project_level_validation_20261007/tables'
effects = source(muscle / 'project_family_effects.csv')
meta = source(muscle / 'project_family_overall_meta.csv')
directions = source(muscle / 'candidate_family_direction.csv')
loo = source(muscle / 'project_family_leave_one_out.csv')


def make(width, height, letter, title):
    fig = plt.figure(figsize=(width * MM, height * MM))
    fig.text(.025, .965, letter, weight='bold', fontsize=11, va='top', color=INK)
    fig.text(.082 if width < 120 else .055, .965, title, fontsize=9, va='top', color=INK)
    return fig


def save(fig, name):
    for ext in ['pdf', 'svg', 'png']:
        fig.savefig(PANELS / f'{name}.{ext}', dpi=450, facecolor='white')
    plt.close(fig)


def clean(ax, grid='x'):
    ax.spines[['top', 'right']].set_visible(False)
    ax.tick_params(length=3, color=INK)
    ax.set_axisbelow(True)
    ax.grid(axis=grid, color='#E9EAEB', linewidth=.5)


def cohort_legend(ax, **kwargs):
    handles = [Line2D([0], [0], color=c, marker='o', lw=1.2, ms=3.5, label=n)
               for _, n, c in COHORTS]
    ax.legend(handles=handles, frameon=False, fontsize=8.5, **kwargs)


def panel1a():
    fig = make(180, 95, 'A', 'From established ADL difficulty to the next observed state')
    ax = fig.add_axes([.025, .045, .95, .84])
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis('off')

    def box(x, y, w, h, text, colour, fill='#FFFFFF', size=8.5):
        ax.add_patch(Rectangle((x, y), w, h, facecolor=fill, edgecolor=colour, lw=.9))
        ax.text(x + w / 2, y + h / 2, text, ha='center', va='center', fontsize=size, linespacing=1.25)

    for i, (cohort, name, colour) in enumerate(COHORTS):
        row = flow[flow.cohort == cohort].iloc[0]
        y = .58 if i == 0 else .10
        box(.00, y, .255, .30,
            ('HRS primary cohort\n2016 Venous Blood Study\nFunction: 2016–2022' if i == 0 else
             'CHARLS\nexploratory comparison\n2011 blood panel\nFunction: 2011, 2013,\n2015, 2018'), colour, size=8.2)
        box(.295, y, .285, .30,
            f'Index visit: ≥1 ADL difficulty\n{int(row.difficulty_start_intervals):,} adjacent-wave intervals\n'
            f'{int(row.unknown_next_wave):,} unknown outcomes excluded', colour)
        box(.625, y, .365, .30,
            f'Eligible analysis risk set\n{int(row.final_participants):,} participants · {int(row.final_intervals):,} intervals\n'
            'Recovery / Persistent difficulty / Death', colour, size=8.4)
        for x0, x1 in [(.255, .295), (.58, .625)]:
            ax.add_patch(FancyArrowPatch((x0, y + .15), (x1, y + .15), arrowstyle='-|>',
                                        mutation_scale=10, lw=.8, color=colour))
        ax.text(.595, y - .025, f'{int(row.excluded_marker_covariate_weight):,} further exclusions: biomarkers, covariates or weights',
                ha='center', va='top', fontsize=7.6, color=INK)
    ax.text(.005, .97, 'Age ≥60 years at biomarker baseline; five harmonised ADL items', fontsize=8.5)
    ax.text(.005, .00, 'Public skeletal muscle transcriptomes: unlinked biological context (Figure 3)', fontsize=8)
    save(fig, 'Figure1A_risk_set_flow')


def panel1b():
    fig = make(90, 83, 'B', 'Observed states at the next wave')
    ax = fig.add_axes([.17, .27, .80, .58])
    states = [('recovery', 'Recovery', RECOVERY),
              ('persistent_difficulty', 'Persistent difficulty', PERSISTENT), ('death', 'Death', DEATH)]
    for x, (cohort, name, _) in enumerate(COHORTS):
        r = flow[flow.cohort == cohort].iloc[0]
        bottom = 0
        for key, label, colour in states:
            value = 100 * r[key] / r.final_intervals
            ax.bar(x, value, bottom=bottom, width=.55, color=colour, edgecolor='white', lw=.5)
            ax.text(x, bottom + value / 2, f'n = {int(r[key]):,}', ha='center', va='center', fontsize=8,
                    color='white' if key == 'death' else INK)
            bottom += value
        ax.text(x, 104, f'n = {int(r.final_intervals):,}', ha='center', fontsize=8)
    ax.set_xticks([0, 1], ['HRS', 'CHARLS'])
    ax.set_ylim(0, 114); ax.set_yticks([0, 25, 50, 75, 100])
    ax.set_ylabel('Eligible intervals (%)')
    ax.set_xlim(-.55, 1.55)
    clean(ax, 'y')
    fig.legend(handles=[Patch(facecolor=c, label=l) for _, l, c in states],
               loc='lower center', bbox_to_anchor=(.50, .00), frameon=False, fontsize=8, ncol=1)
    save(fig, 'Figure1B_observed_states')


def panel1c():
    fig = make(90, 83, 'C', 'Distribution of the signed difference')
    ax = fig.add_axes([.19, .20, .77, .66])
    lo, hi = histogram.left.min(), histogram.right.max()
    for cohort, name, colour in COHORTS:
        data = histogram[histogram.cohort == cohort]
        edges = np.r_[data.left.to_numpy(), data.right.iloc[-1]]
        ax.stairs(data.density, edges, color=colour, linewidth=1.3, label=name)
        ax.stairs(data.density, edges, color=colour, alpha=.08, fill=True)
    ax.set_xlim(lo, hi)
    ytop = ax.get_ylim()[1]
    for i, (cohort, name, colour) in enumerate(COHORTS):
        q = quantiles[quantiles.cohort == cohort].iloc[0]
        y = ytop * (.77 if i == 0 else .52)
        for key in ['P10', 'P50', 'P90']:
            ax.plot([q[key], q[key]], [0, y], color=colour, ls=':', lw=.7)
            ax.text(q[key], y, key, color=colour, fontsize=7.5, ha='center', va='bottom')
    ax.set_xlabel('Signed eGFRcys − eGFRcr difference\n(mL/min/1.73 m²)', fontsize=8)
    ax.set_ylabel('Density')
    ax.legend(frameon=False, fontsize=8, loc='upper left')
    clean(ax, 'y')
    save(fig, 'Figure1C_exposure_distribution')


def panel2a():
    fig = make(90, 88, 'A', 'Recovery and secondary observed states')
    ax = fig.add_axes([.37, .18, .59, .65])
    for outcome, y in [('Recovery', 3), ('Persistent difficulty', 2), ('Death', 1)]:
        for j, (cohort, name, colour) in enumerate(COHORTS):
            r = clinical[(clinical.cohort == cohort) & (clinical.outcome == outcome)].iloc[0]
            yi = y + (.14 if j == 0 else -.14)
            ax.errorbar(r.odds_ratio, yi, xerr=[[r.odds_ratio-r.ci_lower], [r.ci_upper-r.odds_ratio]],
                        fmt='o', color=colour, ms=4, capsize=2, lw=1.2)
    ax.axvline(1, color=DEATH, ls='--', lw=.8)
    ax.set_yticks([3, 2, 1], ['Recovery', 'Persistent\ndifficulty', 'Death'])
    ax.set_ylim(.35, 3.55); ax.set_xscale('log'); ax.set_xlim(.6, 1.35)
    ax.xaxis.set_major_locator(FixedLocator([.6, .8, 1, 1.2]))
    ax.xaxis.set_major_formatter(FixedFormatter(['0.6', '0.8', '1.0', '1.2']))
    ax.xaxis.set_minor_locator(FixedLocator([]))
    ax.set_xlabel('Adjacent-wave OR (95% CI)')
    ax.text(.03, .01, 'Death: I² = 83.3%', transform=ax.transAxes, fontsize=7.8)
    cohort_legend(ax, loc='upper left', bbox_to_anchor=(-.45, 1.18), ncol=2)
    clean(ax)
    save(fig, 'Figure2A_clinical_associations')


def panel2b(letter='B', filename='Figure2B_existing_recovery_splines'):
    fig = make(90, 88, letter, 'Shape of the recovery association')
    ax = fig.add_axes([.19, .20, .77, .65])
    for i, (cohort, name, colour) in enumerate(COHORTS):
        r = curves[curves.cohort == cohort]
        ax.fill_between(r.signed_difference, 100*r.ci_lower, 100*r.ci_upper, color=colour, alpha=.12)
        ax.plot(r.signed_difference, 100*r.probability, color=colour, lw=1.4, label=name)
        q = quantiles[quantiles.cohort == cohort].iloc[0]
        for key in ['P10', 'P90']:
            y = np.interp(q[key], r.signed_difference, 100*r.probability)
            ax.plot(q[key], y, 'o', color=colour, ms=3)
            ax.annotate(key, (q[key], y), xytext=(0, -13 if i == 0 else 9), textcoords='offset points',
                        fontsize=7.5, color=colour, ha='center')
    ax.set_ylim(0, 65)
    ax.set_ylabel('Standardised recovery probability (%)')
    ax.set_xlabel('Signed eGFRcys − eGFRcr difference\n(mL/min/1.73 m²)', fontsize=8)
    ax.legend(frameon=False, loc='upper left', fontsize=8)
    clean(ax, 'y')
    save(fig, filename)


def panel2c():
    fig = make(90, 88, 'C', 'Recovery across the exposure distribution')
    ax = fig.add_axes([.23, .22, .73, .61])
    for i, (cohort, name, colour) in enumerate(COHORTS):
        r = enh[(enh.cohort == cohort) & (enh.analysis == 'standardized_p90_vs_p10') &
                (enh.outcome == 'recovery')].iloc[0]
        y = 1-i
        est, lo, hi = r.estimate*100, r.ci_lower*100, r.ci_upper*100
        ax.errorbar(est, y, xerr=[[est-lo], [hi-est]], fmt='o', color=colour, ms=4.5, capsize=2, lw=1.4)
        label = f'{est:.1f} ({lo:.1f} to {hi:.1f})'.replace('-', '−')
        ax.text(est, y+.25, label, color=colour, fontsize=8, ha='center')
    ax.axvline(0, color=DEATH, ls='--', lw=.8)
    ax.set_yticks([1, 0], ['HRS', 'CHARLS']); ax.set_ylim(-.5, 1.6); ax.set_xlim(-2, 16)
    ax.set_xlabel('P90 − P10 recovery probability\ndifference (percentage points)')
    clean(ax)
    save(fig, 'Figure2C_probability_contrasts')


def panel2d():
    fig = make(90, 88, 'D', 'Recovery in sensitivity analyses')
    ax = fig.add_axes([.48, .18, .48, .66])
    specs = [('primary_mean', 'full', 'Main model\n(arithmetic mean eGFR)'),
             ('primary_combined', 'full', 'Combined creatinine–\ncystatin C eGFR'),
             ('combined_base', 'markers', 'Matched cases:\ncombined eGFR'),
             ('marker_adjusted', 'markers', 'Matched cases:\n+ blood markers')]
    for i, (analysis, sample, label) in enumerate(specs):
        for j, (cohort, name, colour) in enumerate(COHORTS):
            if sample == 'markers':
                target_sample = 'bun_hgb_alb' if cohort.startswith('HRS') else 'bun_hgb'
                target_analysis = ('combined_bun_hgb_alb' if cohort.startswith('HRS') else 'combined_bun_hgb') if analysis == 'marker_adjusted' else analysis
            else:
                target_sample = sample; target_analysis = analysis
            r = renal[(renal.cohort == cohort) & (renal.outcome == 'recovery') &
                      (renal.analysis == target_analysis) & (renal['sample'] == target_sample) &
                      (renal.term == 'delta10')].iloc[0]
            y = 4-i + (.13 if j == 0 else -.13)
            ax.errorbar(r.odds_ratio, y, xerr=[[r.odds_ratio-r.ci_lower], [r.ci_upper-r.odds_ratio]],
                        fmt='o', color=colour, ms=3.8, capsize=2, lw=1.2)
    ax.axvline(1, color=DEATH, ls='--', lw=.8)
    ax.set_yticks([4, 3, 2, 1], [s[2] for s in specs]); ax.tick_params(axis='y', labelsize=7.8)
    ax.set_ylim(.45, 4.55); ax.set_xlim(.97, 1.28)
    ax.set_xticks([1, 1.1, 1.2], ['1.0', '1.1', '1.2'])
    ax.set_xlabel('Recovery OR (95% CI)')
    cohort_legend(ax, loc='upper left', bbox_to_anchor=(-.82, 1.17), ncol=2)
    clean(ax)
    save(fig, 'Figure2D_recovery_sensitivity')


def family_name(name):
    return 'GSE111017 family' if '111017' in name else name


def panel3a():
    fig = make(180, 68, 'A', 'OXPHOS scores in adverse muscle phenotypes')
    ax = fig.add_axes([.26, .25, .71, .60])
    for i, r in effects.iterrows():
        y = 4-i
        ax.errorbar(r.beta, y, xerr=[[r.beta-r.ci_low], [r.ci_high-r.beta]], fmt='o',
                    color=LOWER, ms=4, capsize=2, lw=1.3)
    r = meta.iloc[0]
    ax.errorbar(r.beta, 0, xerr=[[r.beta-r.ci_low], [r.ci_high-r.beta]], fmt='s',
                color=HRS, ms=4.5, capsize=2, lw=1.5)
    ax.axvline(0, color=DEATH, ls='--', lw=.8)
    ax.set_yticks([4, 3, 2, 1, 0], [family_name(n) for n in effects.study_family] + ['Pooled (4 families)'])
    ax.set_ylim(-.5, 4.5); ax.set_xlim(-1.5, .7)
    ax.set_xlabel('Standardised OXPHOS score contrast, beta (95% CI)')
    fig.text(.97, .018, 'Pooled beta −0.298 (95% CI −0.512 to −0.084); I² = 0%',
             ha='right', fontsize=8)
    clean(ax)
    save(fig, 'Figure3A_OXPHOS_family_forest')


def panel3b():
    fig = make(90, 86, 'B', 'Candidate-gene direction by family')
    ax = fig.add_axes([.20, .28, .66, .55])
    families = effects.study_family.tolist()
    genes = ['COX5A', 'MDH2', 'SDHA', 'GOT2']
    for i, gene in enumerate(genes):
        for j, family in enumerate(families):
            r = directions[(directions.symbol == gene) & (directions.study_family == family)].iloc[0]
            colour = LOWER if r.mean_logFC < 0 else HIGHER
            ax.scatter(j, 3-i, s=135, marker='o' if i < 3 else 'D', color=colour, edgecolor='white', lw=.7)
        ax.text(3.7, 3-i, '4/4' if i < 3 else '3/4', va='center', fontsize=8)
    ax.set_xticks(range(4), ['GSE111017\nfamily', 'GSE117525', 'GSE144304', 'GSE174106'],
                  rotation=45, ha='right', fontsize=7.7)
    ax.set_yticks([3, 2, 1, 0], genes)
    ax.set_xlim(-.5, 4.25); ax.set_ylim(-.5, 3.6)
    ax.text(3.8, 3.6, 'Lower\nfamilies', ha='center', fontsize=7.3)
    ax.tick_params(length=0)
    for spine in ax.spines.values(): spine.set_visible(False)
    fig.legend(handles=[Line2D([0], [0], marker='o', color='none', markerfacecolor=LOWER, label='Lower'),
                        Line2D([0], [0], marker='o', color='none', markerfacecolor=HIGHER, label='Higher')],
               loc='lower center', bbox_to_anchor=(.48, .06), frameon=False, ncol=2, fontsize=8)
    fig.legend(handles=[Line2D([0], [0], marker='o', color='none', markerfacecolor=INK, label='4/4 lower'),
                        Line2D([0], [0], marker='D', color='none', markerfacecolor=INK, label='3/4 lower')],
               loc='lower center', bbox_to_anchor=(.48, .005), frameon=False, ncol=2, fontsize=7.8)
    save(fig, 'Figure3B_candidate_directions')


def panel3c():
    fig = make(90, 86, 'C', 'Leave-one-family-out sensitivity')
    ax = fig.add_axes([.47, .20, .49, .64])
    for i, r in loo.iterrows():
        ax.errorbar(r.beta, 3-i, xerr=[[r.beta-r.ci_low], [r.ci_high-r.beta]],
                    fmt='o', color=LOWER, ms=4, capsize=2, lw=1.2)
    ax.axvline(0, color=DEATH, ls='--', lw=.8)
    ax.set_yticks([3, 2, 1, 0], ['Omit ' + family_name(n).replace(' family', '') for n in loo.omitted_project_family])
    ax.tick_params(axis='y', labelsize=8)
    ax.set_ylim(-.5, 3.6); ax.set_xlim(-.8, .2); ax.set_xticks([-.8, -.4, 0])
    ax.set_xlabel('Pooled beta (95% CI)')
    clean(ax)
    save(fig, 'Figure3C_leave_one_family_out')


def assemble(stem, width, height, placements, folder=None):
    folder = folder or V1 / 'figures'
    folder.mkdir(parents=True, exist_ok=True)
    doc = fitz.open()
    page = doc.new_page(width=width * 72 / 25.4, height=height * 72 / 25.4)
    for name, x, y, w, h in placements:
        panel_pdf = fitz.open(PANELS / f'{name}.pdf')
        rect = fitz.Rect(x*72/25.4, y*72/25.4, (x+w)*72/25.4, (y+h)*72/25.4)
        page.show_pdf_page(rect, panel_pdf, 0)
        panel_pdf.close()
    doc.save(folder / f'{stem}.pdf', deflate=True)
    page.get_pixmap(dpi=600, alpha=False).save(folder / f'{stem}.png')
    with Image.open(folder / f'{stem}.png') as im:
        im.save(folder / f'{stem}.tif', compression='tiff_lzw', dpi=(600, 600))
    doc.close()


def supplementary_shape():
    panel2b('A', 'FigureS1A_existing_recovery_splines')
    fig = make(180, 85, 'B', 'Exploratory interaction tests for recovery')
    ax = fig.add_axes([.32, .18, .64, .65])
    names = ['age', 'sex', 'BMI', 'CKD', 'baseline_severity']
    labels = ['Age', 'Sex', 'Body mass index', 'Arithmetic mean eGFR <60', 'Index ADL difficulty count']
    # Preserve the original 30-test multiplicity family, even though this panel
    # displays its recovery subset only.
    allp = interactions.p.to_numpy()
    order = np.argsort(allp)
    adjusted = np.empty(len(allp))
    adjusted[order] = np.minimum(1, np.minimum.accumulate((allp[order] * len(allp) / np.arange(1, len(allp)+1))[::-1])[::-1])
    z = interactions.copy(); z['p_bh_display'] = adjusted
    for j, (cohort, name, colour) in enumerate(COHORTS):
        rows = z[(z.cohort == cohort) & (z.outcome == 'recovery')].set_index('modifier').loc[names]
        y = np.arange(5)[::-1] + (.12 if j == 0 else -.12)
        ax.plot(rows.p_bh_display, y, 'o', color=colour, ms=4, label=name)
    ax.axvline(.05, color=DEATH, ls='--', lw=.8)
    ax.set_yticks(np.arange(5)[::-1], labels)
    ax.set_xlim(0, 1); ax.set_xlabel(r'BH-adjusted interaction $P$ (30-test family)')
    ax.legend(frameon=False, fontsize=8, ncol=2, loc='upper right')
    clean(ax)
    save(fig, 'FigureS1_recovery_interactions')
    assemble('FigureS1_recovery_shape_and_interactions', 180, 173,
             [('FigureS1A_existing_recovery_splines', 45, 0, 90, 88),
              ('FigureS1_recovery_interactions', 0, 88, 180, 85)], V1 / 'supplementary_figures')


def supplementary_nonlinearity():
    fig = make(180, 90, 'A', 'Nonlinearity tests for all observed next-wave states')
    ax = fig.add_axes([.12, .18, .84, .65])
    for j, (cohort, name, colour) in enumerate(COHORTS):
        r = nonlinear[nonlinear.cohort == cohort].set_index('outcome').loc[['recovery', 'persistent_difficulty', 'death']]
        xs = np.arange(3) + (-.08 if j == 0 else .08)
        ax.plot(xs, -np.log10(r.p_nonlinear), 'o-', color=colour, ms=4, lw=1, label=name)
        for x, p in zip(xs, r.p_nonlinear):
            label = '$P$ <0.001' if p < .001 else f'$P$ = {p:.3f}'
            ax.annotate(label, (x, -np.log10(p)), xytext=(0, 8 if j == 0 else -14),
                        textcoords='offset points', fontsize=8, color=colour, ha='center')
    ax.set_xticks([0, 1, 2], ['Recovery', 'Persistent difficulty', 'Death'])
    ax.set_xlim(-.4, 2.4); ax.set_ylim(-.13, max(-np.log10(nonlinear.p_nonlinear))+.45)
    ax.set_ylabel(r'$-\log_{10}$($P$ for nonlinearity)')
    ax.legend(frameon=False, fontsize=8, loc='upper left', ncol=2)
    clean(ax, 'y')
    save(fig, 'FigureS3_nonlinearity_tests')
    assemble('FigureS3_nonlinearity_tests', 180, 90,
             [('FigureS3_nonlinearity_tests', 0, 0, 180, 90)], V1 / 'supplementary_figures')


