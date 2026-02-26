"""
Plot Figure 2: Treatment outcomes under standard of care.

Shows overtreatment vs undertreatment disaggregated by:
    - Diagnostic pathway (GUD syndromic, ANC screening, secondary rash, tertiary)
    - Sex
    - HIV status
"""

import sciris as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as pl
from matplotlib.gridspec import GridSpec
from utils import set_font

RESULTS_DIR = 'results'
FIGURES_DIR = 'figures'

# Pathway display names and colors
PATHWAY_LABELS = {
    'gud_syndromic': 'GUD\nsyndromic',
    'anc_screen': 'ANC\nscreening',
    'secondary_rash': 'Secondary\nrash',
    'tertiary': 'Tertiary',
}
PATHWAY_COLORS = {
    'gud_syndromic': '#e41a1c',
    'anc_screen': '#377eb8',
    'secondary_rash': '#ff7f00',
    'tertiary': '#984ea3',
}

# Outcome colors
OC_COLORS = {
    'success': '#4daf4a',       # Green - correctly treated
    'unnecessary': '#e41a1c',   # Red - overtreated
    'failure': '#ff7f00',       # Orange - treatment failed
    'missed': '#984ea3',        # Purple - missed diagnoses
}

# Sex colors
F_COLOR = '#d46e9c'
M_COLOR = '#4a90d9'

# HIV colors
HIV_POS_COLOR = '#e41a1c'
HIV_NEG_COLOR = '#377eb8'


def load_data(scenario='soc'):
    """Load treatment outcomes data"""
    df = sc.loadobj(f'{RESULTS_DIR}/treatment_outcomes_{scenario}.df')
    return df


def pivot_annual(df, start_year=2000, end_year=2035):
    """Get annual averages across parameter sets for each metric"""
    df = df[(df.year >= start_year) & (df.year <= end_year)].copy()
    df['year_int'] = df.year.astype(int)
    annual = df.groupby(['scenario', 'year_int', 'metric']).value.agg(['mean', 'std', 'count']).reset_index()
    return annual


def get_metric(df, metric, start_year=2000, end_year=2035):
    """Get annual mean for a specific metric"""
    sub = df[(df.metric == metric) & (df.year_int >= start_year) & (df.year_int <= end_year)]
    return sub.set_index('year_int')['mean']


def plot_treatment_time_series(df, ax, start_year=2000, end_year=2035):
    """Panel A: Treatment outcomes over time (stacked area)"""
    pathways = ['gud_syndromic', 'anc_screen', 'secondary_rash', 'tertiary']

    # Total treated = success + unnecessary + failure across all pathways
    for i, pw in enumerate(pathways):
        success = get_metric(df, f'{pw}_success', start_year, end_year)
        unnecessary = get_metric(df, f'{pw}_unnecessary', start_year, end_year)
        total = success + unnecessary
        years = total.index

        ax.plot(years, total, color=PATHWAY_COLORS[pw], linewidth=2,
                label=PATHWAY_LABELS[pw].replace('\n', ' '))

    ax.set_ylabel('Treatments per year')
    ax.set_title('Treatments by pathway')
    ax.legend(frameon=False, fontsize=11, loc='upper left')
    ax.set_xlim(start_year, end_year)
    ax.set_ylim(bottom=0)
    sc.SIticks(ax, axis='y')


def plot_overtreatment_by_pathway(df, ax, start_year=2000, end_year=2035):
    """Panel B: Stacked bars showing treated vs overtreated by pathway"""
    pathways = ['gud_syndromic', 'anc_screen', 'secondary_rash', 'tertiary']
    x = np.arange(len(pathways))
    width = 0.35

    success_vals = []
    unnecessary_vals = []
    for pw in pathways:
        s = get_metric(df, f'{pw}_success', start_year, end_year)
        u = get_metric(df, f'{pw}_unnecessary', start_year, end_year)
        success_vals.append(s.mean())
        unnecessary_vals.append(u.mean())

    bars_s = ax.bar(x - width/2, success_vals, width, label='Correctly treated',
                     color=OC_COLORS['success'], alpha=0.85)
    bars_u = ax.bar(x + width/2, unnecessary_vals, width, label='Overtreated',
                     color=OC_COLORS['unnecessary'], alpha=0.85)

    ax.set_xticks(x)
    ax.set_xticklabels([PATHWAY_LABELS[pw] for pw in pathways])
    ax.set_ylabel('Annual average')
    ax.set_title('Overtreatment by pathway')
    ax.legend(frameon=False, fontsize=11)
    ax.set_ylim(bottom=0)
    sc.SIticks(ax, axis='y')


def plot_missed_by_pathway(df, ax, start_year=2000, end_year=2035):
    """Panel C: Missed diagnoses (false negatives) by pathway"""
    pathways = ['gud_syndromic', 'anc_screen', 'secondary_rash']  # No testing for tertiary
    x = np.arange(len(pathways))

    missed_vals = []
    for pw in pathways:
        m = get_metric(df, f'{pw}_missed', start_year, end_year)
        missed_vals.append(m.mean())

    bars = ax.bar(x, missed_vals, color=[PATHWAY_COLORS[pw] for pw in pathways], alpha=0.85)
    for bar, val in zip(bars, missed_vals):
        if val > 0:
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
                    f'{val:.0f}', ha='center', va='bottom', fontsize=12)

    ax.set_xticks(x)
    ax.set_xticklabels([PATHWAY_LABELS[pw] for pw in pathways])
    ax.set_ylabel('Annual average')
    ax.set_title('Missed diagnoses\n(false negatives)')
    ax.set_ylim(bottom=0)
    sc.SIticks(ax, axis='y')


def plot_outcomes_by_sex(df, ax, start_year=2000, end_year=2035):
    """Panel D: Treatment outcomes by sex"""
    pathways = ['gud_syndromic', 'anc_screen', 'secondary_rash', 'tertiary']
    outcomes = ['success', 'unnecessary']
    x = np.arange(len(pathways))
    width = 0.2

    for i, (sex, color, label) in enumerate([(('_f', F_COLOR, 'Female')), ('_m', M_COLOR, 'Male')]):
        success = []
        unnecessary = []
        for pw in pathways:
            s = get_metric(df, f'{pw}_success{sex}', start_year, end_year)
            u = get_metric(df, f'{pw}_unnecessary{sex}', start_year, end_year)
            success.append(s.mean())
            unnecessary.append(u.mean())

        offset = (i - 0.5) * width * 2
        ax.bar(x + offset - width/2, success, width, color=color, alpha=0.85,
               label=f'{label} correct' if i == 0 else f'{label} correct')
        ax.bar(x + offset + width/2, unnecessary, width, color=color, alpha=0.35,
               label=f'{label} overtreated' if i == 0 else f'{label} overtreated',
               edgecolor=color, linewidth=1)

    ax.set_xticks(x)
    ax.set_xticklabels([PATHWAY_LABELS[pw] for pw in pathways])
    ax.set_ylabel('Annual average')
    ax.set_title('Treatment outcomes by sex')
    ax.legend(frameon=False, fontsize=10, ncol=2)
    ax.set_ylim(bottom=0)
    sc.SIticks(ax, axis='y')


def plot_outcomes_by_hiv(df, ax, start_year=2000, end_year=2035):
    """Panel E: Treatment outcomes by HIV status"""
    pathways = ['gud_syndromic', 'anc_screen', 'secondary_rash', 'tertiary']
    x = np.arange(len(pathways))
    width = 0.2

    for i, (suffix, color, label) in enumerate([('_hivpos', HIV_POS_COLOR, 'HIV+'), ('_hivneg', HIV_NEG_COLOR, 'HIV\u2212')]):
        success = []
        unnecessary = []
        for pw in pathways:
            s = get_metric(df, f'{pw}_success{suffix}', start_year, end_year)
            u = get_metric(df, f'{pw}_unnecessary{suffix}', start_year, end_year)
            success.append(s.mean())
            unnecessary.append(u.mean())

        offset = (i - 0.5) * width * 2
        ax.bar(x + offset - width/2, success, width, color=color, alpha=0.85,
               label=f'{label} correct')
        ax.bar(x + offset + width/2, unnecessary, width, color=color, alpha=0.35,
               label=f'{label} overtreated', edgecolor=color, linewidth=1)

    ax.set_xticks(x)
    ax.set_xticklabels([PATHWAY_LABELS[pw] for pw in pathways])
    ax.set_ylabel('Annual average')
    ax.set_title('Treatment outcomes by HIV status')
    ax.legend(frameon=False, fontsize=10, ncol=2)
    ax.set_ylim(bottom=0)
    sc.SIticks(ax, axis='y')


def plot_overtreatment_rate(df, ax, start_year=2000, end_year=2035):
    """Panel F: Overtreatment rate (% of treatments that are unnecessary) by pathway over time"""
    pathways = ['gud_syndromic', 'anc_screen', 'secondary_rash']

    for pw in pathways:
        success = get_metric(df, f'{pw}_success', start_year, end_year)
        unnecessary = get_metric(df, f'{pw}_unnecessary', start_year, end_year)
        total = success + unnecessary
        rate = unnecessary / total.where(total > 0) * 100
        years = rate.index
        ax.plot(years, rate, color=PATHWAY_COLORS[pw], linewidth=2,
                label=PATHWAY_LABELS[pw].replace('\n', ' '))

    ax.set_ylabel('Overtreatment rate (%)')
    ax.set_title('Overtreatment rate by pathway')
    ax.legend(frameon=False, fontsize=11)
    ax.set_xlim(start_year, end_year)
    ax.set_ylim(0, 100)


if __name__ == '__main__':

    scenario = 'soc'
    df = load_data(scenario)
    annual = pivot_annual(df, start_year=2000, end_year=2035)

    set_font(size=18)
    fig = pl.figure(figsize=(22, 14))
    gs = GridSpec(2, 3, left=0.06, right=0.98, bottom=0.06, top=0.93,
                  wspace=0.28, hspace=0.38)

    ax = fig.add_subplot(gs[0, 0])
    plot_treatment_time_series(annual, ax)

    ax = fig.add_subplot(gs[0, 1])
    plot_overtreatment_by_pathway(annual, ax)

    ax = fig.add_subplot(gs[0, 2])
    plot_missed_by_pathway(annual, ax)

    ax = fig.add_subplot(gs[1, 0])
    plot_outcomes_by_sex(annual, ax)

    ax = fig.add_subplot(gs[1, 1])
    plot_outcomes_by_hiv(annual, ax)

    ax = fig.add_subplot(gs[1, 2])
    plot_overtreatment_rate(annual, ax)

    for i, ax in enumerate(fig.axes):
        ax.text(-0.10, 1.06, chr(65 + i), transform=ax.transAxes,
                fontsize=28, fontweight='bold', va='top')

    pl.savefig(f'{FIGURES_DIR}/fig2_treatment_outcomes_soc.png', dpi=200, bbox_inches='tight')
    print(f'Saved {FIGURES_DIR}/fig2_treatment_outcomes_soc.png')
