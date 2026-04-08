"""
Slide panel 1: "The Problem" — unnecessary syphilis treatments are large and growing.

Single figure (3.83 × 3.83 in) — grouped bars: historical (2015–2025) vs projected
(2026–2035), split by symptomatic (GUD syndromic) vs asymptomatic (ANC+KP+PLHIV screening).
"""

import sciris as sc
import numpy as np
import matplotlib.pyplot as pl
import matplotlib.patches as mpatches
import sys; sys.path.insert(0, '..')
from utils import set_font

RESULTS_DIR = '../results'
FIGURES_DIR = '.'

SYMP_PATHWAYS  = ['gud_syndromic']
ASYMP_PATHWAYS = ['anc_screen', 'kp_screen', 'plhiv_screen']

GROUPS = [
    ('Symptomatic\n(syndromic\nmanagement)', SYMP_PATHWAYS,  '#e41a1c'),
    ('Asymptomatic\n(treponemal\nscreening)',  ASYMP_PATHWAYS, '#377eb8'),
]

HIST = (2015, 2025)
PROJ = (2026, 2035)


def load_data():
    return sc.loadobj(f'{RESULTS_DIR}/treatment_outcomes_soc.df').copy()


def sum_pathways(df, pathways, start, end):
    sub = df[(df['year'] >= start) & (df['year'] <= end)]
    correct = sum(sub[f'{pw}_success'].mean() for pw in pathways)
    unnecessary = sum(sub[f'{pw}_unnecessary'].mean() for pw in pathways)
    return correct, unnecessary


def plot():
    set_font(size=12)
    fig, ax = pl.subplots(figsize=(3.83, 3.83), facecolor='white')
    fig.subplots_adjust(left=0.22, right=0.95, top=0.95, bottom=0.22)

    df = load_data()

    n = len(GROUPS)
    x = np.arange(n)
    width = 0.33
    gap = 0.04

    for i, (label, pathways, color) in enumerate(GROUPS):
        hc, hu = sum_pathways(df, pathways, *HIST)
        pc, pu = sum_pathways(df, pathways, *PROJ)

        # Historical bar (left)
        xh = x[i] - width/2 - gap
        ax.bar(xh, hc, width, color=color, alpha=0.90)
        ax.bar(xh, hu, width, bottom=hc, color=color, alpha=0.25,
               hatch='////', edgecolor=color, linewidth=0.5)
        tot_h = hc + hu
        ot_h = hu / tot_h * 100 if tot_h > 0 else 0

        # Projected bar (right)
        xp = x[i] + width/2 + gap
        ax.bar(xp, pc, width, color=color, alpha=0.90)
        ax.bar(xp, pu, width, bottom=pc, color=color, alpha=0.25,
               hatch='////', edgecolor=color, linewidth=0.5)
        tot_p = pc + pu
        ot_p = pu / tot_p * 100 if tot_p > 0 else 0

        # Labels above bars
        label_offset = 800
        ax.text(xh, tot_h + label_offset, f'{tot_h/1e3:.0f}K\n({ot_h:.0f}% unn.)',
                ha='center', va='bottom', fontsize=8.5, color='#555555',
                fontweight='bold', linespacing=1.1)
        ax.text(xp, tot_p + label_offset, f'{tot_p/1e3:.0f}K\n({ot_p:.0f}% unn.)',
                ha='center', va='bottom', fontsize=8.5, color='#555555',
                fontweight='bold', linespacing=1.1)

        # Period labels below bars
        ax.text(xh, -2200, f'{HIST[0]}–\n{HIST[1]}', ha='center', va='top',
                fontsize=7, color='#999999', linespacing=1.1)
        ax.text(xp, -2200, f'{PROJ[0]}–\n{PROJ[1]}', ha='center', va='top',
                fontsize=7, color='#999999', linespacing=1.1)

    ax.set_xticks(x)
    ax.set_xticklabels([g[0] for g in GROUPS], fontsize=9)
    ax.set_ylabel('Treatments per year')
    ax.set_ylim(bottom=0)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(axis='x', length=0, pad=28)
    sc.SIticks(ax, axis='y')

    legend_handles = [
        mpatches.Patch(facecolor='#666666', alpha=0.90, label='Correctly treated'),
        mpatches.Patch(facecolor='#666666', alpha=0.25, hatch='////',
                       edgecolor='#666666', label='Unnecessary'),
    ]
    ax.legend(handles=legend_handles, loc='upper center', fontsize=8.5,
              frameon=False, handlelength=1.2, handletextpad=0.5, ncol=1)

    pl.savefig(f'{FIGURES_DIR}/slide_problem.png', dpi=300, facecolor='white', bbox_inches=None)
    print(f'Saved {FIGURES_DIR}/slide_problem.png')
    pl.close(fig)


if __name__ == '__main__':
    plot()
