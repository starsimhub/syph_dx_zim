"""
Slide panel 2: "The Solution" — waterfall showing how new diagnostics remove
unnecessary treatments, broken out by product and use case.

Single figure (3.83 × 3.83 in). Starts with SOC total, steps down by each
diagnostic use case, ends at residual.
"""

import sciris as sc
import numpy as np
import matplotlib.pyplot as pl
import matplotlib.patches as mpatches
import sys; sys.path.insert(0, '..')
from utils import set_font

RESULTS_DIR = '../results'
FIGURES_DIR = '.'

SYMP  = ['gud_syndromic']
ASYMP = ['anc_screen', 'kp_screen', 'plhiv_screen']
ALL   = SYMP + ASYMP
PERIOD = (2026, 2035)


def load_scenario(scenario):
    df = sc.loadobj(f'{RESULTS_DIR}/treatment_outcomes_{scenario}.df').copy()
    sub = df[(df['year'] >= PERIOD[0]) & (df['year'] <= PERIOD[1])]
    return sub


def get_totals(sub, pathways):
    c = sum(sub[f'{pw}_success'].mean() for pw in pathways)
    u = sum(sub[f'{pw}_unnecessary'].mean() for pw in pathways)
    return c, u


def plot():
    set_font(size=12)
    fig, ax = pl.subplots(figsize=(3.83, 3.83), facecolor='white')
    fig.subplots_adjust(left=0.20, right=0.95, top=0.95, bottom=0.14)

    soc_sub  = load_scenario('soc')
    gud_sub  = load_scenario('gud')
    both_sub = load_scenario('both')

    # Totals
    soc_c, soc_u = get_totals(soc_sub, ALL)
    _, soc_u_symp = get_totals(soc_sub, SYMP)
    _, gud_u_symp = get_totals(gud_sub, SYMP)
    both_c, both_u = get_totals(both_sub, ALL)

    # Per-pathway NT drops
    nt_drops = {}
    for pw in ASYMP:
        soc_pw_u = soc_sub[f'{pw}_unnecessary'].mean()
        both_pw_u = both_sub[f'{pw}_unnecessary'].mean()
        nt_drops[pw] = soc_pw_u - both_pw_u

    # Waterfall steps
    correct = soc_c
    drop_gud = soc_u_symp - gud_u_symp
    drop_anc = nt_drops['anc_screen']
    drop_kp = nt_drops['kp_screen']
    drop_plhiv = nt_drops['plhiv_screen']

    # Running total of unnecessary (ordered by drop size)
    levels = [
        soc_u,
        soc_u - drop_gud,
        soc_u - drop_gud - drop_plhiv,
        soc_u - drop_gud - drop_plhiv - drop_kp,
        soc_u - drop_gud - drop_plhiv - drop_kp - drop_anc,
    ]

    drops = [None, drop_gud, drop_plhiv, drop_kp, drop_anc]
    labels = ['SOC', 'GUD\nPOC test', 'PLHIV\nPOC NT', 'KP\nPOC NT', 'ANC\nPOC NT']

    x = np.arange(len(labels))
    width = 0.58

    # Colors
    col_unn = '#E8453C'
    col_correct = '#1B9E77'
    col_drop = '#3B82F6'

    # Correct treatment baseline across all bars
    for xi in x:
        ax.bar(xi, correct, width, color=col_correct, alpha=0.90, zorder=2)

    # SOC bar: full unnecessary stack
    ax.bar(x[0], levels[0], width, bottom=correct, color=col_unn, alpha=0.85, zorder=2)

    # Drop bars (hanging)
    for i in range(1, len(x)):
        ax.bar(x[i], drops[i], width, bottom=correct + levels[i],
               color=col_drop, alpha=0.80, zorder=2)

    # Remaining unnecessary in final bar
    ax.bar(x[-1], levels[-1], width, bottom=correct, color=col_unn, alpha=0.85, zorder=2)

    # Connector lines
    for i in range(len(x) - 1):
        y_connect = correct + levels[i + 1]
        ax.plot([x[i] + width/2, x[i+1] - width/2], [y_connect, y_connect],
                color='#999999', linewidth=0.8, linestyle='--', zorder=1)

    # Drop labels inside bars
    for i in range(1, len(x)):
        mid_y = correct + levels[i] + drops[i] / 2
        if drops[i] > 5000:
            ax.text(x[i], mid_y, f'–{drops[i]/1e3:.0f}K',
                    ha='center', va='center', fontsize=9, fontweight='bold',
                    color='white', zorder=3)
        else:
            # Small drop — label above the bar
            ax.text(x[i], correct + levels[i-1] + 800, f'–{drops[i]/1e3:.1f}K',
                    ha='center', va='bottom', fontsize=7.5, fontweight='bold',
                    color=col_drop, zorder=3)

    # SOC total label above bar
    ax.text(x[0], correct + levels[0] + 600, f'{(correct + levels[0])/1e3:.0f}K',
            ha='center', va='bottom', fontsize=9, fontweight='bold', color='#333333')

    # SOC unnecessary label inside red area
    ax.text(x[0], correct + levels[0] / 2, f'{levels[0]/1e3:.0f}K',
            ha='center', va='center', fontsize=9, fontweight='bold',
            color='white', zorder=3)

    # Remaining label inside the red area
    red_mid = correct + levels[-1] / 2
    ax.text(x[-1], red_mid, f'{levels[-1]/1e3:.0f}K',
            ha='center', va='center', fontsize=9, fontweight='bold',
            color='white', zorder=3)

    # Correct treatment label in first bar
    ax.text(x[0], correct / 2, f'{correct/1e3:.0f}K\ncorrectly\ntreated',
            ha='center', va='center', fontsize=7.5, color='white', fontweight='bold',
            linespacing=1.15, zorder=3)

    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylabel('Treatments per year\n(average over 2026–2035)')
    ax.set_ylim(bottom=0)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(axis='x', length=0)
    sc.SIticks(ax, axis='y')

    legend_handles = [
        mpatches.Patch(facecolor=col_correct, alpha=0.90, label='Correctly treated'),
        mpatches.Patch(facecolor=col_unn, alpha=0.85, label='Unnecessary'),
        mpatches.Patch(facecolor=col_drop, alpha=0.80, label='Avoided by new dx'),
    ]
    ax.legend(handles=legend_handles, loc='upper right', fontsize=7.5,
              frameon=False, handlelength=1.2, handletextpad=0.4)

    pl.savefig(f'{FIGURES_DIR}/slide_solution.png', dpi=300, facecolor='white', bbox_inches=None)
    print(f'Saved {FIGURES_DIR}/slide_solution.png')
    pl.close(fig)


if __name__ == '__main__':
    plot()
