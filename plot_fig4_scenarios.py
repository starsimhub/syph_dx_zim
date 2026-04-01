"""
Plot Figure 4: Scenario comparison — impact of novel diagnostics.

Panels:
    A: Total treatments per year by scenario (adult pathways)
    B: Adult overtreatment rate (%) over time
    C: Average overtreatment rate (2028–2039) by scenario
"""

import sciris as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as pl
from matplotlib.gridspec import GridSpec
from utils import set_font, get_metric

RESULTS_DIR = 'results'
FIGURES_DIR = 'figures'

SCENARIOS = ['soc', 'gud', 'conf', 'both']

SCENARIO_LABELS = {
    'soc': 'Standard\nof care',
    'gud': 'GUD\nPOC test',
    'conf': 'POC NT active\ninfection diagnostic',
    'both': 'Both',
}
SCENARIO_COLORS = {
    'soc': '#555555',
    'gud': '#e41a1c',
    'conf': '#377eb8',
    'both': '#4daf4a',
}

ADULT_PATHWAYS = ['gud_syndromic', 'anc_screen', 'kp_screen', 'plhiv_screen']


def load_all_scenarios(scenarios=None):
    if scenarios is None:
        scenarios = SCENARIOS
    dfs = {}
    for scen in scenarios:
        try:
            dfs[scen] = sc.loadobj(f'{RESULTS_DIR}/treatment_outcomes_{scen}.df')
        except FileNotFoundError:
            print(f'WARNING: {scen} results not found, skipping')
    return dfs



def get_pathway_total(ann, metric_suffix, pathways, start_year=2025, end_year=2040):
    """Sum a metric across specified pathways"""
    total = None
    for pw in pathways:
        s = get_metric(ann, f'{pw}_{metric_suffix}', start_year, end_year)
        if total is None:
            total = s
        else:
            total = total.add(s, fill_value=0)
    return total


def plot_treatments_ts(dfs, ax, start_year=2025, end_year=2039):
    """Panel A: Adult treatments per year by scenario"""
    for scen in SCENARIOS:
        if scen not in dfs:
            continue
        adult = get_pathway_total(dfs[scen], 'treated', ADULT_PATHWAYS, start_year, end_year)
        color = SCENARIO_COLORS[scen]
        label = SCENARIO_LABELS[scen].replace('\n', ' ')
        ax.plot(adult.index, adult, color=color, linewidth=2.5, label=label)

    ax.set_ylabel('Treatments per year')
    ax.set_title('(A) Syphilis treatments\nby scenario')
    ax.legend(frameon=False, fontsize=14)
    ax.set_xlim(start_year, end_year)
    ax.set_ylim(bottom=0)
    sc.SIticks(ax, axis='y')


def plot_overtreatment_ts(dfs, ax, start_year=2025, end_year=2039):
    """Panel B: Adult overtreatment rate (%) over time"""
    for scen in SCENARIOS:
        if scen not in dfs:
            continue
        treated = get_pathway_total(dfs[scen], 'treated', ADULT_PATHWAYS, start_year, end_year)
        unnecessary = get_pathway_total(dfs[scen], 'unnecessary', ADULT_PATHWAYS, start_year, end_year)
        rate = (unnecessary / treated * 100).clip(upper=100)

        color = SCENARIO_COLORS[scen]
        label = SCENARIO_LABELS[scen].replace('\n', ' ')
        ax.plot(rate.index, rate, color=color, linewidth=2.5, label=label)

    ax.set_ylabel('Overtreatment rate (%)')
    ax.set_title('(B) Adult overtreatment rate\nby scenario')
    ax.legend(frameon=False, fontsize=14)
    ax.set_xlim(start_year, end_year)
    ax.set_ylim(0, 100)


def plot_overtreatment_bars(dfs, ax, start_year=2028, end_year=2039):
    """Panel C: Average overtreatment rate (2028–2039) by scenario"""
    bar_labels = [SCENARIO_LABELS[s].replace('\n', '\n') for s in SCENARIOS]
    bar_colors = [SCENARIO_COLORS[s] for s in SCENARIOS]

    rates = []
    for scen in SCENARIOS:
        if scen not in dfs:
            rates.append(0)
            continue
        treated = get_pathway_total(dfs[scen], 'treated', ADULT_PATHWAYS, start_year, end_year)
        unnecessary = get_pathway_total(dfs[scen], 'unnecessary', ADULT_PATHWAYS, start_year, end_year)
        total_treated = treated.sum()
        total_unnecessary = unnecessary.sum()
        rates.append(total_unnecessary / total_treated * 100 if total_treated > 0 else 0)

    x = np.arange(len(SCENARIOS))
    bars = ax.bar(x, rates, color=bar_colors, alpha=0.85, edgecolor='white', linewidth=0.5, width=0.6)
    for bar, rate in zip(bars, rates):
        ax.text(bar.get_x() + bar.get_width()/2, max(bar.get_height(), 0) + 1,
                f'{rate:.0f}%', ha='center', va='bottom', fontsize=14, fontweight='bold')

    ax.set_xticks(x)
    ax.set_xticklabels(bar_labels)
    ax.set_ylabel('Overtreatment rate (%)')
    ax.set_title(f'(C) Overtreatment rate\n{start_year}\u2013{end_year}')
    ax.set_ylim(0, 110)


def print_summary_table(dfs, start_year=2028, end_year=2039):
    """Print a summary table of key metrics across scenarios"""
    print(f'\n{"="*90}')
    print(f'Summary: {start_year}-{end_year} annual averages')
    print(f'{"="*90}')
    header = f'{"Scenario":<18} {"Treated/yr":>12} {"Correct/yr":>12} {"Overtreat/yr":>12} {"OT rate":>8} {"Pathways":>14}'
    print(header)
    print('-' * len(header))

    for scen in SCENARIOS:
        if scen not in dfs:
            continue
        pws = ADULT_PATHWAYS
        pw_label = '(adult)'
        treated = get_pathway_total(dfs[scen], 'treated', pws, start_year, end_year).mean()
        success = get_pathway_total(dfs[scen], 'success', pws, start_year, end_year).mean()
        unnecessary = get_pathway_total(dfs[scen], 'unnecessary', pws, start_year, end_year).mean()
        ot_rate = unnecessary / treated * 100 if treated > 0 else 0
        label = SCENARIO_LABELS[scen].replace('\n', ' ')
        print(f'{label:<18} {treated:>12,.0f} {success:>12,.0f} {unnecessary:>12,.0f} {ot_rate:>7.1f}% {pw_label:>14}')

    print(f'{"="*90}\n')


if __name__ == '__main__':

    dfs = load_all_scenarios()
    if not dfs:
        print('No scenario results found. Run run_scenarios.py first.')
        exit()

    print_summary_table(dfs)

    set_font(size=20)
    fig = pl.figure(figsize=(22, 8))
    gs = GridSpec(1, 3, left=0.05, right=0.98, bottom=0.10, top=0.88,
                  wspace=0.18)

    ax = fig.add_subplot(gs[0, 0])
    plot_treatments_ts(dfs, ax)

    ax = fig.add_subplot(gs[0, 1])
    plot_overtreatment_ts(dfs, ax)

    ax = fig.add_subplot(gs[0, 2])
    plot_overtreatment_bars(dfs, ax)

    pl.savefig(f'{FIGURES_DIR}/fig4_scenario_comparison.png', dpi=200, bbox_inches='tight')
    print(f'Saved {FIGURES_DIR}/fig4_scenario_comparison.png')
