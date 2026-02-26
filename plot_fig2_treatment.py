"""
Plot Figure 2: Treatment outcomes under standard of care.

Panels:
    A: Treatments by pathway over time
    B: Stacked bars — correct + overtreated + missed by pathway
    C: Overtreatment rate (%) by pathway
    D: Treatment outcomes by sex
"""

import sciris as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as pl
from matplotlib.gridspec import GridSpec
from utils import set_font

RESULTS_DIR = 'results'
FIGURES_DIR = 'figures'

PATHWAY_LABELS = {
    'gud_syndromic': 'GUD\nsyndromic',
    'anc_screen': 'ANC\nscreening',
    'secondary_rash': 'Secondary\nrash',
    'kp_screen': 'KP dual\nRDT',
}
PATHWAY_COLORS = {
    'gud_syndromic': '#e41a1c',
    'anc_screen': '#377eb8',
    'secondary_rash': '#ff7f00',
    'kp_screen': '#984ea3',
}

OC_COLORS = {
    'success': '#4daf4a',
    'unnecessary': '#e41a1c',
    'missed': '#984ea3',
}

F_COLOR = '#d46e9c'
M_COLOR = '#4a90d9'

PATHWAYS = ['gud_syndromic', 'anc_screen', 'secondary_rash', 'kp_screen']


def load_data(scenario='soc'):
    return sc.loadobj(f'{RESULTS_DIR}/treatment_outcomes_{scenario}.df')


def pivot_annual(df, start_year=2000, end_year=2035):
    df = df[(df.year >= start_year) & (df.year <= end_year)].copy()
    df['year_int'] = df.year.astype(int)
    return df.groupby(['scenario', 'year_int', 'metric']).value.agg(['mean', 'std', 'count']).reset_index()


def get_metric(df, metric, start_year=2000, end_year=2035):
    sub = df[(df.metric == metric) & (df.year_int >= start_year) & (df.year_int <= end_year)]
    return sub.set_index('year_int')['mean']


def plot_treatment_time_series(df, ax, start_year=2000, end_year=2035):
    """Panel A: Treatments by pathway over time"""
    for pw in PATHWAYS:
        success = get_metric(df, f'{pw}_success', start_year, end_year)
        unnecessary = get_metric(df, f'{pw}_unnecessary', start_year, end_year)
        total = success + unnecessary
        ax.plot(total.index, total, color=PATHWAY_COLORS[pw], linewidth=2,
                label=PATHWAY_LABELS[pw].replace('\n', ' '))

    ax.set_ylabel('Treatments per year')
    ax.set_title('Treatments by pathway')
    ax.legend(frameon=False, fontsize=13, loc='upper left')
    ax.set_xlim(start_year, end_year)
    ax.set_ylim(bottom=0)
    sc.SIticks(ax, axis='y')


def plot_stacked_outcomes(df, ax, start_year=2000, end_year=2035):
    """Panel B: Stacked bars — correct (green) + overtreated (red) + missed (purple) by pathway"""
    x = np.arange(len(PATHWAYS))
    width = 0.6

    success_vals, unnecessary_vals, missed_vals = [], [], []
    for pw in PATHWAYS:
        success_vals.append(get_metric(df, f'{pw}_success', start_year, end_year).mean())
        unnecessary_vals.append(get_metric(df, f'{pw}_unnecessary', start_year, end_year).mean())
        missed_vals.append(get_metric(df, f'{pw}_missed', start_year, end_year).mean())

    s = np.array(success_vals)
    u = np.array(unnecessary_vals)
    m = np.array(missed_vals)

    ax.bar(x, s, width, label='Correctly treated', color=OC_COLORS['success'], alpha=0.85)
    ax.bar(x, u, width, bottom=s, label='Overtreated', color=OC_COLORS['unnecessary'], alpha=0.85)
    ax.bar(x, m, width, bottom=s + u, label='Missed', color=OC_COLORS['missed'], alpha=0.85)

    # Total labels
    for i in range(len(PATHWAYS)):
        total = s[i] + u[i] + m[i]
        if total > 0:
            ax.text(x[i], total + total * 0.02, f'{total:,.0f}',
                    ha='center', va='bottom', fontsize=12)

    ax.set_xticks(x)
    ax.set_xticklabels([PATHWAY_LABELS[pw] for pw in PATHWAYS])
    ax.set_ylabel('Annual average')
    ax.set_title('Treatment outcomes by pathway')
    ax.legend(frameon=False, fontsize=13)
    ax.set_ylim(bottom=0)
    sc.SIticks(ax, axis='y')


def plot_overtreatment_rate_bars(df, ax, start_year=2000, end_year=2035):
    """Panel C: Overtreatment rate (%) by pathway"""
    x = np.arange(len(PATHWAYS))

    rates = []
    for pw in PATHWAYS:
        s = get_metric(df, f'{pw}_success', start_year, end_year).mean()
        u = get_metric(df, f'{pw}_unnecessary', start_year, end_year).mean()
        total = s + u
        rates.append(u / total * 100 if total > 0 else 0)

    bars = ax.bar(x, rates, color=[PATHWAY_COLORS[pw] for pw in PATHWAYS], alpha=0.85)
    for bar, rate in zip(bars, rates):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 1,
                f'{rate:.0f}%', ha='center', va='bottom', fontsize=16, fontweight='bold')

    ax.set_xticks(x)
    ax.set_xticklabels([PATHWAY_LABELS[pw] for pw in PATHWAYS])
    ax.set_ylabel('Overtreatment rate (%)')
    ax.set_title('Overtreatment rate by pathway')
    ax.set_ylim(0, 110)


def plot_outcomes_by_sex(df, ax, start_year=2000, end_year=2035):
    """Panel D: Treatment outcomes by sex"""
    x = np.arange(len(PATHWAYS))
    width = 0.2

    for i, (sex, color, label) in enumerate([('_f', F_COLOR, 'Female'), ('_m', M_COLOR, 'Male')]):
        success, unnecessary = [], []
        for pw in PATHWAYS:
            success.append(get_metric(df, f'{pw}_success{sex}', start_year, end_year).mean())
            unnecessary.append(get_metric(df, f'{pw}_unnecessary{sex}', start_year, end_year).mean())

        offset = (i - 0.5) * width * 2
        ax.bar(x + offset - width/2, success, width, color=color, alpha=0.85,
               label=f'{label} correct')
        ax.bar(x + offset + width/2, unnecessary, width, color=color, alpha=0.35,
               label=f'{label} overtreated', edgecolor=color, linewidth=1)

    ax.set_xticks(x)
    ax.set_xticklabels([PATHWAY_LABELS[pw] for pw in PATHWAYS])
    ax.set_ylabel('Annual average')
    ax.set_title('Treatment outcomes by sex')
    ax.legend(frameon=False, fontsize=12, ncol=2)
    ax.set_ylim(bottom=0)
    sc.SIticks(ax, axis='y')


if __name__ == '__main__':

    scenario = 'soc'
    df = load_data(scenario)
    annual = pivot_annual(df, start_year=2000, end_year=2035)

    set_font(size=18)
    fig = pl.figure(figsize=(18, 14))
    gs = GridSpec(2, 2, left=0.07, right=0.98, bottom=0.06, top=0.93,
                  wspace=0.28, hspace=0.38)

    ax = fig.add_subplot(gs[0, 0])
    plot_treatment_time_series(annual, ax)

    ax = fig.add_subplot(gs[0, 1])
    plot_stacked_outcomes(annual, ax)

    ax = fig.add_subplot(gs[1, 0])
    plot_overtreatment_rate_bars(annual, ax)

    ax = fig.add_subplot(gs[1, 1])
    plot_outcomes_by_sex(annual, ax)

    for i, ax in enumerate(fig.axes):
        ax.text(-0.10, 1.06, chr(65 + i), transform=ax.transAxes,
                fontsize=28, fontweight='bold', va='top')

    pl.savefig(f'{FIGURES_DIR}/fig2_treatment_outcomes_soc.png', dpi=200, bbox_inches='tight')
    print(f'Saved {FIGURES_DIR}/fig2_treatment_outcomes_soc.png')
