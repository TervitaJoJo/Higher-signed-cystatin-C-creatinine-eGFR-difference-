from pathlib import Path
import hashlib
import importlib.util
import json
import shutil
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FixedLocator, FixedFormatter
from matplotlib.patches import Rectangle
import fitz
from PIL import Image

import os
ROOT = Path(os.environ["REPRO_WORK"])
AUDIT = ROOT / "analysis/audit/bmc_clinical_focus_20261011"
AUDIT.mkdir(parents=True, exist_ok=True)
import plot_base as draw
DATA = draw.renal
source_rows = []


def model(cohort, analysis, sample):
    rows = DATA[(DATA.cohort == cohort) & (DATA.outcome == 'recovery') &
                (DATA.term == 'delta10') & (DATA.analysis == analysis) & (DATA['sample'] == sample)]
    assert len(rows) == 1
    row = rows.iloc[0]
    source_rows.append(row.to_dict())
    return row


def panel2a():
    fig = draw.make(180, 72, 'A', 'Recovery and secondary observed states')
    ax = fig.add_axes([.22, .23, .74, .57])
    for outcome, y in [('Recovery', 3), ('Persistent difficulty', 2), ('Death', 1)]:
        for j, (cohort, name, colour) in enumerate(COHORTS):
            row = draw.clinical[(draw.clinical.cohort == cohort) & (draw.clinical.outcome == outcome)].iloc[0]
            ax.errorbar(row.odds_ratio, y + (.14 if j == 0 else -.14),
                        xerr=[[row.odds_ratio-row.ci_lower], [row.ci_upper-row.odds_ratio]],
                        fmt='o' if j == 0 else 's', color=colour, ms=4, capsize=2, lw=1.3)
    ax.axvline(1, color=draw.DEATH, ls='--', lw=.8)
    ax.set_yticks([3, 2, 1], ['Recovery', 'Persistent difficulty', 'Death'])
    ax.set_ylim(.45, 3.55)
    ax.set_xscale('log')
    ax.set_xlim(.6, 1.35)
    ax.xaxis.set_major_locator(FixedLocator([.6, .8, 1, 1.2]))
    ax.xaxis.set_major_formatter(FixedFormatter(['0.6', '0.8', '1.0', '1.2']))
    ax.xaxis.set_minor_locator(FixedLocator([]))
    ax.set_xlabel('Adjacent-wave OR (95% CI)')
    handles = [Line2D([0], [0], color=c, marker='o' if i == 0 else 's', lw=0, ms=4,
                      label='HRS primary cohort' if i == 0 else 'CHARLS exploratory comparison')
               for i, (_, _, c) in enumerate(COHORTS)]
    fig.legend(handles=handles, loc='upper center', bbox_to_anchor=(.63, .925),
               frameon=False, ncol=2, fontsize=8.5)
    fig.text(.025, .025, 'Death: I² = 83.3%; no stable pooled interpretation', fontsize=8)
    draw.clean(ax)
    draw.save(fig, 'Figure2A_clinical_associations')


def panel3a():
    fig = draw.make(180, 88, 'A', 'Alternative filtration adjustment')
    for index, (cohort, name, colour) in enumerate(COHORTS):
        ax = fig.add_axes([.12 + index*.47, .24, .34, .56])
        rows = [model(cohort, analysis, 'full') for analysis in ['primary_mean', 'primary_combined']]
        estimate = np.array([r.odds_ratio for r in rows])
        lower = np.array([r.ci_lower for r in rows])
        upper = np.array([r.ci_upper for r in rows])
        ax.plot([0, 1], estimate, color=colour, lw=1.3, zorder=3)
        for x, row in enumerate(rows):
            ax.errorbar(x, estimate[x], yerr=[[estimate[x]-lower[x]], [upper[x]-estimate[x]]],
                        fmt='o' if index == 0 else 's', color=colour, ms=5,
                        markerfacecolor='white' if x == 0 else colour, capsize=4, lw=1.3, zorder=4)
            ax.annotate(f'{row.odds_ratio:.3f}', (x, estimate[x]), xytext=(9, 5),
                        textcoords='offset points', fontsize=9, color=colour)
        ax.axhline(1, color=draw.DEATH, ls='--', lw=.8)
        ax.set_ylim(.97, 1.28)
        ax.set_xlim(-.35, 1.42)
        ax.set_yticks([1, 1.1, 1.2], ['1.0', '1.1', '1.2'])
        ax.set_xticks([0, 1], ['Arithmetic\nmean eGFR', 'Combined\neGFR'])
        ax.tick_params(axis='x', labelsize=8.5)
        if index == 0:
            ax.set_ylabel('Recovery OR (95% CI)')
        draw.clean(ax, 'y')
        fig.text(.29 + index*.47, .865,
                 'HRS primary cohort' if index == 0 else 'CHARLS exploratory comparison',
                 ha='center', fontsize=9, color=colour)
        fig.text(.29 + index*.47, .055,
                 f'{int(rows[0].n_people):,} participants; {int(rows[0].n_intervals):,} intervals',
                 ha='center', fontsize=8)
    draw.save(fig, 'Figure3A_filtration_adjustment_paths')


def p_label(p):
    return r'$\mathit{P}$ <0.001' if p < .001 else r'$\mathit{P}$' + f' = {p:.3f}'


def panel3b():
    fig = draw.make(180, 90, 'B', 'Blood-marker adjustment within matched complete cases')
    ax = fig.add_axes([.025, .065, .95, .79])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis('off')
    left = .28
    middle = .64
    ax.text((left+middle)/2, .99, 'Combined eGFR\n+ clinical covariates',
            ha='center', va='top', fontsize=9)
    ax.text((middle+1)/2, .99, '+ blood markers\n(same complete cases)',
            ha='center', va='top', fontsize=9)
    for y in [.79, .405, .025]:
        ax.plot([0, 1], [y, y], color='#C6C9CB', lw=.7)
    ax.plot([left, left], [.025, .79], color='#E2E4E5', lw=.6)
    ax.plot([middle, middle], [.025, .79], color='#E2E4E5', lw=.6)
    for index, (cohort, name, colour) in enumerate(COHORTS):
        sample = 'bun_hgb_alb' if index == 0 else 'bun_hgb'
        added = 'combined_bun_hgb_alb' if index == 0 else 'combined_bun_hgb'
        rows = [model(cohort, 'combined_base', sample), model(cohort, added, sample)]
        assert all(rows[0][key] == rows[1][key] for key in ['n_people', 'n_intervals', 'events'])
        top = .755-index*.385
        ax.text(.01, top-.025, name, ha='left', va='top', fontsize=10, color=colour, fontweight='bold')
        ax.text(.01, top-.135,
                f'{int(rows[0].n_people):,} participants\n{int(rows[0].n_intervals):,} intervals\n'
                f'{int(rows[0].events):,} recoveries', ha='left', va='top', fontsize=8.5, linespacing=1.45)
        for j, row in enumerate(rows):
            x = (left+middle)/2 if j == 0 else (middle+1)/2
            ax.text(x, top-.01, f'OR {row.odds_ratio:.3f}', ha='center', va='top',
                    fontsize=11, color=colour, fontweight='bold')
            ax.text(x, top-.13, f'95% CI {row.ci_lower:.3f}–{row.ci_upper:.3f}',
                    ha='center', va='top', fontsize=9)
            ax.text(x, top-.25, p_label(row.p), ha='center', va='top', fontsize=9)
    draw.save(fig, 'Figure3B_matched_case_estimate_matrix')


def assemble_svg(stem, width, height, placements):
    from lxml import etree
    namespace = 'http://www.w3.org/2000/svg'
    svg = etree.Element('{%s}svg' % namespace, nsmap={None: namespace},
                        width=f'{width}mm', height=f'{height}mm', viewBox=f'0 0 {width} {height}')
    for index, (name, x, y, w, h) in enumerate(placements):
        child = etree.parse(str(draw.PANELS / f'{name}.svg')).getroot()
        replacements = {}
        for element in child.iter():
            if element.get('id'):
                old = element.get('id')
                replacements[old] = f'p{index}_{old}'
                element.set('id', replacements[old])
        for element in child.iter():
            for key, value in list(element.attrib.items()):
                for old, new in replacements.items():
                    value = value.replace(f'url(#{old})', f'url(#{new})')
                    if value == f'#{old}': value = f'#{new}'
                element.set(key, value)
        for key, value in {'x': x, 'y': y, 'width': w, 'height': h}.items():
            child.set(key, str(value))
        svg.append(child)
    etree.ElementTree(svg).write(str(V1 / 'figures' / f'{stem}.svg'), encoding='utf-8', xml_declaration=True)


