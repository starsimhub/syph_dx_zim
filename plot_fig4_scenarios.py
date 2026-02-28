"""
Plot Figure 3: Scenario comparison — impact of novel diagnostics.

Panels:
    A: Total treatments per year — adult (solid) and newborn (dashed) lines
    B: Adult overtreatment rate (%) over time (excl. CS — adult diagnostics only)
    C: Overtreatment rate bars — adults for GUD/conf/both, newborn for CS
"""

import sciris as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as pl
from matplotlib.gridspec import GridSpec
from utils import set_font

RESULTS_DIR = 'results'
FIGURES_DIR = 'figures'

SCENARIOS = ['soc', 'gud', 'conf', 'both', 'cs']
ADULT_SCENARIOS = ['soc', 'gud', 'conf', 'both']  # Scenarios that affect adult pathways

SCENARIO_LABELS = {
    'soc': 'Standard\nof care',
    'gud': 'GUD\nPOC test',
    'conf': 'Confirmatory\ntest',
    'both': 'Both',
    'cs': 'Newborn\nPOC CS',
}
SCENARIO_COLORS = {
    'soc': '#555555',
    'gud': '#e41a1c',
    'conf': '#377eb8',
    'both': '#4daf4a',
    'cs': '#ff7f00',
}

ADULT_PATHWAYS = ['gud_syndromic', 'anc_screen', 'kp_screen', 'plhiv_screen']
NEWBORN_PATHWAYS = ['newborn']
ALL_PATHWAYS = ADULT_PATHWAYS + NEWBORN_PATHWAYS


def load_all_scenarios():
    dfs = {}
    for scen in SCENARIOS:
        try:
            dfs[scen] = sc.loadobj(f'{RESULTS_DIR}/treatment_outcomes_{scen}.df')
        except FileNotFoundError:
            print(f'WARNING: {scen} results not found, skipping')
    return dfs


def pivot_annual(df, start_year=2025, end_year=2040):
    df = df[(df.year >= start_year) & (df.year <= end_year)].copy()
    df['year_int'] = df.year.astype(int)
    return df.groupby(['scenario', 'year_int', 'metric']).value.agg(['mean', 'std']).reset_index()


def get_metric(ann, metric, start_year=2025, end_year=2040):
    sub = ann[(ann.metric == metric) & (ann.year_int >= start_year) & (ann.year_int <= end_year)]
    return sub.set_index('year_int')['mean']


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
    """Panel A: Adult treatments (solid) and newborn treatments (dashed) by scenario"""
    for scen in SCENARIOS:
        if scen not in dfs:
            continue
        ann = pivot_annual(dfs[scen], start_year, end_year)
        adult = get_pathway_total(ann, 'treated', ADULT_PATHWAYS, start_year, end_year)
        newborn = get_pathway_total(ann, 'treated', NEWBORN_PATHWAYS, start_year, end_year)

        color = SCENARIO_COLORS[scen]
        label = SCENARIO_LABELS[scen].replace('\n', ' ')

        # Only plot adult lines for adult scenarios, newborn for cs
        if scen in ADULT_SCENARIOS:
            ax.plot(adult.index, adult, color=color, linewidth=2.5, label=label)
        if scen in ['soc', 'cs']:
            ls = '-' if scen == 'cs' else '--'
            lbl = label if scen == 'cs' else f'{label} (newborn)'
            ax.plot(newborn.index, newborn, color=SCENARIO_COLORS[scen], linewidth=2,
                    linestyle=ls, label=lbl, alpha=0.8)

    ax.set_ylabel('Treatments per year')
    ax.set_title('(A) Syphilis treatments\nby scenario')
    ax.legend(frameon=False, fontsize=14)
    ax.set_xlim(start_year, end_year)
    ax.set_ylim(bottom=0)
    sc.SIticks(ax, axis='y')


def plot_overtreatment_ts(dfs, ax, start_year=2025, end_year=2039):
    """Panel B: Adult overtreatment rate (%) over time — adult scenarios only"""
    for scen in ADULT_SCENARIOS:
        if scen not in dfs:
            continue
        ann = pivot_annual(dfs[scen], start_year, end_year)
        treated = get_pathway_total(ann, 'treated', ADULT_PATHWAYS, start_year, end_year)
        unnecessary = get_pathway_total(ann, 'unnecessary', ADULT_PATHWAYS, start_year, end_year)
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
    """Panel C: Overtreatment rate bars — adult pathways for adult scenarios, newborn for CS"""
    x = np.arange(len(SCENARIOS))
    rates = []

    for scen in SCENARIOS:
        if scen not in dfs:
            rates.append(0)
            continue
        ann = pivot_annual(dfs[scen], start_year, end_year)
        if scen == 'cs':
            # CS: compute newborn overtreatment rate
            treated = get_pathway_total(ann, 'treated', NEWBORN_PATHWAYS, start_year, end_year)
            unnecessary = get_pathway_total(ann, 'unnecessary', NEWBORN_PATHWAYS, start_year, end_year)
        else:
            # Adult scenarios: compute adult overtreatment rate
            treated = get_pathway_total(ann, 'treated', ADULT_PATHWAYS, start_year, end_year)
            unnecessary = get_pathway_total(ann, 'unnecessary', ADULT_PATHWAYS, start_year, end_year)

        total_treated = treated.sum()
        total_unnecessary = unnecessary.sum()
        rate = total_unnecessary / total_treated * 100 if total_treated > 0 else 0
        rates.append(rate)

    colors = [SCENARIO_COLORS[s] for s in SCENARIOS]
    bars = ax.bar(x, rates, color=colors, alpha=0.85, edgecolor='white', linewidth=0.5, width=0.6)
    for bar, rate in zip(bars, rates):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 1,
                f'{rate:.0f}%', ha='center', va='bottom', fontsize=14, fontweight='bold')

    ax.set_xticks(x)
    ax.set_xticklabels([SCENARIO_LABELS[s] for s in SCENARIOS])
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
        ann = pivot_annual(dfs[scen], start_year, end_year)
        if scen == 'cs':
            pws = NEWBORN_PATHWAYS
            pw_label = '(newborn)'
        else:
            pws = ADULT_PATHWAYS
            pw_label = '(adult)'
        treated = get_pathway_total(ann, 'treated', pws, start_year, end_year).mean()
        success = get_pathway_total(ann, 'success', pws, start_year, end_year).mean()
        unnecessary = get_pathway_total(ann, 'unnecessary', pws, start_year, end_year).mean()
        ot_rate = unnecessary / treated * 100 if treated > 0 else 0
        label = SCENARIO_LABELS[scen].replace('\n', ' ')
        print(f'{label:<18} {treated:>12,.0f} {success:>12,.0f} {unnecessary:>12,.0f} {ot_rate:>7.1f}% {pw_label:>14}')

    print(f'{"="*90}\n')


def plot_congenital_cascade_comparison(dfs, cs, ax_soc, ax_cs, start_year=2028, end_year=2039, anc_attendance=0.90):
    """Row 2: Side-by-side congenital cascades for SOC vs CS scenario"""
    cols = cs.columns.get_level_values(0)

    # MTC transmissions (same denominator for both — ANC pathway unchanged)
    mtc_total = sum(
        cs.loc[(cs.index >= start_year) & (cs.index <= end_year), (f'transmission_by_stage.new_mtc_{s}', '50%')].mean()
        for s in ['primary', 'secondary', 'early', 'late', 'tertiary']
        if f'transmission_by_stage.new_mtc_{s}' in cols
    )

    for scen_key, ax, title_prefix in [('soc', ax_soc, '(D) SOC'), ('cs', ax_cs, '(E) Newborn POC')]:
        if scen_key not in dfs:
            continue
        df = dfs[scen_key]
        sub = df[(df.year >= start_year) & (df.year <= end_year)]

        anc_success = sub[sub.metric == 'anc_screen_success'].groupby('year').value.mean().mean()
        anc_missed = sub[sub.metric == 'anc_screen_missed'].groupby('year').value.mean().mean()
        anc_failure = sub[sub.metric == 'anc_screen_failure'].groupby('year').value.mean().mean()
        nb_success = sub[sub.metric == 'newborn_success'].groupby('year').value.mean().mean()

        mothers_screened_active = anc_success + anc_missed + anc_failure
        fetal_prevented = anc_success * 0.90
        total_at_risk = mtc_total + fetal_prevented

        f = 100 / total_at_risk if total_at_risk > 0 else 1
        attend_anc = 100 * anc_attendance
        screened = mothers_screened_active * f
        detected = (anc_success + anc_failure) * f
        treated = anc_success * f
        fetus_cured = fetal_prevented * f
        remaining = 100 - fetus_cured
        nb_saved = nb_success / mtc_total * remaining if mtc_total > 0 else 0

        steps = [100, attend_anc, screened, detected, treated, fetus_cured]
        labels = ['At-risk\npregnancies', 'Attend\nANC', 'Screened for\nsyphilis',
                  'Test\npositive', 'Mother\ntreated', 'Fetus cured\nin utero']
        loss_labels = ['', 'No ANC', 'Not screened', 'False negative', 'Not treated', 'Treatment failure']

        bar_color = '#377eb8'
        loss_color = '#dddddd'
        y = np.arange(len(steps) + 3)[::-1]

        ax.barh(y[0:len(steps)], steps, color=bar_color, alpha=0.85,
                edgecolor='white', linewidth=0.5, height=0.7)
        for i in range(1, len(steps)):
            ax.barh(y[i], steps[i-1] - steps[i], left=steps[i], color=loss_color, alpha=0.5, height=0.7)

        nb_steps = [remaining, nb_saved]
        nb_labels_list = ['Still at risk\nat birth', 'Newborn\ntreated']
        nb_y = y[len(steps)+1:len(steps)+3]
        ax.barh(nb_y, nb_steps, color='#e41a1c', alpha=0.7, edgecolor='white', linewidth=0.5, height=0.7)
        ax.barh(nb_y[0], 100 - remaining, left=remaining, color=bar_color, alpha=0.3, height=0.7)

        all_steps = list(steps) + [None] + list(nb_steps)
        all_labels = list(labels) + [''] + list(nb_labels_list)
        for i, (step, _) in enumerate(zip(all_steps, all_labels)):
            if step is None:
                continue
            if step >= 1:
                ax.text(step + 1, y[i], f'{step:.0f}', ha='left', va='center', fontsize=16, fontweight='bold')
            elif step > 0:
                ax.text(max(step, 0) + 1, y[i], f'{step:.1f}', ha='left', va='center', fontsize=16, fontweight='bold')

        for i in range(1, len(steps)):
            lost = steps[i-1] - steps[i]
            if lost > 1:
                ax.text(steps[i-1] - 0.5, y[i] + 0.35, f'\u2212{lost:.0f}',
                        ha='right', va='bottom', fontsize=12, color='#888888', style='italic')

        sep_y = (y[len(steps)-1] + y[len(steps)]) / 2
        ax.axhline(y=sep_y, color='grey', linewidth=0.8, linestyle='--', alpha=0.4)

        all_y_labels = list(labels) + [''] + list(nb_labels_list)
        ax.set_yticks(y)
        ax.set_yticklabels(all_y_labels)
        ax.set_xlim(0, 115)
        ax.set_title(f'{title_prefix}: congenital\nprevention cascade')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)


if __name__ == '__main__':

    dfs = load_all_scenarios()
    if not dfs:
        print('No scenario results found. Run run_scenarios.py first.')
        exit()

    print_summary_table(dfs)

    set_font(size=20)
    fig = pl.figure(figsize=(22, 8))
    gs = GridSpec(1, 3, left=0.05, right=0.98, bottom=0.10, top=0.88,
                  wspace=0.25)

    ax = fig.add_subplot(gs[0, 0])
    plot_treatments_ts(dfs, ax)

    ax = fig.add_subplot(gs[0, 1])
    plot_overtreatment_ts(dfs, ax)

    ax = fig.add_subplot(gs[0, 2])
    plot_overtreatment_bars(dfs, ax)

    pl.savefig(f'{FIGURES_DIR}/fig4_scenario_comparison.png', dpi=200, bbox_inches='tight')
    print(f'Saved {FIGURES_DIR}/fig4_scenario_comparison.png')
