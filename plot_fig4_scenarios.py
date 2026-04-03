"""
Plot Figure 4: Scenario comparison — impact of novel diagnostics.

Panels:
    A: Per-pathway SOC vs single-use-case stacked bars (correctly treated +
       overtreated), annotated with % OT change
    B: All 15 non-SOC scenarios ranked by unnecessary treatments avoided vs SOC
    C: Unnecessary treatments over time — SOC, each use case in isolation,
       all confirmatory channels, both products together
    D: Average overtreatment rate (2027–2040) for the same scenarios as C
"""

import sciris as sc
import numpy as np
import matplotlib.pyplot as pl
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Patch
from utils import set_font, get_metric

RESULTS_DIR = 'results'
FIGURES_DIR = 'figures'

# Year ranges
TS_START  = 2025   # time-series window
TS_END    = 2040
BAR_START = 2027   # summary average window
BAR_END   = 2040

ADULT_PATHWAYS = ['gud_syndromic', 'anc_screen', 'kp_screen', 'plhiv_screen']

# Scenarios for panels C and D
TS_SCENARIOS = ['soc', 'gud', 'anc', 'kp', 'plhiv', 'both']

# Non-SOC scenarios for panel B ranking
ALL_NON_SOC = ['gud', 'anc', 'kp', 'plhiv', 'both']

# Pathway → single use-case scenario (for panel A)
PATHWAY_SCENARIOS = {
    'gud_syndromic': 'gud',
    'anc_screen':    'anc',
    'kp_screen':     'kp',
    'plhiv_screen':  'plhiv',
}

PATHWAY_LABELS = {
    'gud_syndromic': 'GUD\nsyndromic',
    'anc_screen':    'ANC\nscreening',
    'kp_screen':     'KP dual\nRDT',
    'plhiv_screen':  'PLHIV\nscreening',
}

SCENARIO_LABELS = {
    'soc':   'SOC',
    'gud':   'GUD\ntest',
    'anc':   'ANC\nconf.',
    'kp':    'KP\nconf.',
    'plhiv': 'PLHIV\nconf.',
    'both':  'All\nproducts',
}

OC_COLORS = {'success': '#4daf4a', 'unnecessary': '#e41a1c'}

TS_COLORS = {
    'soc':   '#555555',
    'gud':   '#e41a1c',
    'anc':   '#377eb8',
    'kp':    '#984ea3',
    'plhiv': '#ff7f00',
    'both':  '#4daf4a',
}


# ── Data helpers ──────────────────────────────────────────────────────────────

def load_scenarios(scenarios):
    dfs = {}
    for scen in scenarios:
        fname = f'{RESULTS_DIR}/treatment_outcomes_{scen}.df'
        try:
            dfs[scen] = sc.loadobj(fname).copy()
        except FileNotFoundError:
            print(f'WARNING: {fname} not found, skipping')
    return dfs


def pathway_metric(df, pathway, suffix, start_year, end_year):
    """Annual mean of {pathway}_{suffix} over [start_year, end_year]."""
    return get_metric(df, f'{pathway}_{suffix}', start_year, end_year).mean()


def total_unnecessary(df, start_year, end_year):
    """Annual time series of unnecessary treatments summed across adult pathways."""
    total = None
    for pw in ADULT_PATHWAYS:
        s = get_metric(df, f'{pw}_unnecessary', start_year, end_year)
        total = s if total is None else total.add(s, fill_value=0)
    return total


def total_treated(df, start_year, end_year):
    """Annual time series of all treatments (success + unnecessary) across adult pathways."""
    total = None
    for pw in ADULT_PATHWAYS:
        s = get_metric(df, f'{pw}_success', start_year, end_year)
        u = get_metric(df, f'{pw}_unnecessary', start_year, end_year)
        t = s.add(u, fill_value=0)
        total = t if total is None else total.add(t, fill_value=0)
    return total


# ── Panel A: per-pathway comparison ──────────────────────────────────────────

def plot_pathway_comparison(dfs, ax, start_year=BAR_START, end_year=BAR_END):
    """SOC vs single use-case stacked bars for each adult pathway."""
    bar_w   = 0.35
    gap     = 0.08   # gap between the two bars in a group
    spacing = 0.25   # extra space between groups

    group_w = 2 * bar_w + gap + spacing
    group_x = np.arange(len(ADULT_PATHWAYS)) * group_w

    for i, pw in enumerate(ADULT_PATHWAYS):
        scen = PATHWAY_SCENARIOS[pw]
        x_soc  = group_x[i]
        x_scen = group_x[i] + bar_w + gap

        for x, key in [(x_soc, 'soc'), (x_scen, scen)]:
            if key not in dfs:
                continue
            s = pathway_metric(dfs[key], pw, 'success',     start_year, end_year)
            u = pathway_metric(dfs[key], pw, 'unnecessary', start_year, end_year)
            ax.bar(x, s, bar_w, color=OC_COLORS['success'],     alpha=0.85)
            ax.bar(x, u, bar_w, color=OC_COLORS['unnecessary'], alpha=0.85, bottom=s)
            total = s + u
            if total > 0:
                ax.text(x, total * 1.02, f'{total:,.0f}',
                        ha='center', va='bottom', fontsize=10)

        # % OT change annotation between the two bars
        if 'soc' in dfs and scen in dfs:
            u_soc  = pathway_metric(dfs['soc'],  pw, 'unnecessary', start_year, end_year)
            s_soc  = pathway_metric(dfs['soc'],  pw, 'success',     start_year, end_year)
            u_scen = pathway_metric(dfs[scen],   pw, 'unnecessary', start_year, end_year)
            if u_soc > 0:
                pct = (u_scen - u_soc) / u_soc * 100
                ax.text((x_soc + x_scen) / 2, (s_soc + u_soc) * 1.18,
                        f'{pct:+.0f}% OT',
                        ha='center', va='bottom', fontsize=10,
                        color='#e41a1c', fontweight='bold')

    tick_x = group_x + bar_w / 2 + gap / 2
    ax.set_xticks(tick_x)
    ax.set_xticklabels([PATHWAY_LABELS[pw] for pw in ADULT_PATHWAYS])
    ax.set_ylabel(f'Annual average ({start_year}–{end_year})')
    ax.set_title('(A) Per-pathway impact\n(left = SOC,  right = use case alone)')
    ax.set_ylim(bottom=0)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    sc.SIticks(ax, axis='y')
    ax.legend(handles=[
        Patch(facecolor=OC_COLORS['success'],     alpha=0.85, label='Correctly treated'),
        Patch(facecolor=OC_COLORS['unnecessary'], alpha=0.85, label='Overtreated'),
    ], frameon=False, fontsize=12, loc='upper right')


# ── Panel B: ranking (two subpanels) ─────────────────────────────────────────

def _ranking_subpanel(ax, soc_u, scenario_list, dfs, title, start_year, end_year,
                      shared_xlim=None, show_ylabel=True):
    """Draw one horizontal-bar ranking subpanel; returns max x value."""
    rows = []
    for scen in scenario_list:
        if scen not in dfs:
            continue
        u       = total_unnecessary(dfs[scen], start_year, end_year).mean()
        avoided = soc_u - u
        color   = TS_COLORS.get(scen, '#888888')
        label   = SCENARIO_LABELS.get(scen, scen).replace('\n', ' ')
        rows.append((label, avoided, color))

    rows.sort(key=lambda r: r[1], reverse=True)

    y    = np.arange(len(rows))
    bars = ax.barh(y, [r[1] for r in rows],
                   color=[r[2] for r in rows], alpha=0.85, height=0.65, edgecolor='white')
    for bar, (_, val, _) in zip(bars, rows):
        ax.text(val + soc_u * 0.005, bar.get_y() + bar.get_height() / 2,
                f'{val:,.0f}', va='center', fontsize=10, fontweight='bold')

    ax.set_yticks(y)
    ax.set_yticklabels([r[0] for r in rows] if show_ylabel else [''] * len(rows), fontsize=11)
    ax.invert_yaxis()
    ax.set_title(title, fontsize=13)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    if shared_xlim is not None:
        ax.set_xlim(0, shared_xlim)
    return max((r[1] for r in rows), default=0)


def plot_ranking(dfs, ax, start_year=BAR_START, end_year=BAR_END):
    """Panel B: ranked horizontal bars for the 4 non-SOC scenarios."""
    if 'soc' not in dfs:
        ax.text(0.5, 0.5, 'SOC data not available', ha='center', transform=ax.transAxes)
        return

    soc_u = total_unnecessary(dfs['soc'], start_year, end_year).mean()
    _ranking_subpanel(ax, soc_u, ALL_NON_SOC, dfs,
                      '(B) Unnecessary treatments avoided vs SOC',
                      start_year, end_year, show_ylabel=True)
    ax.set_xlabel(f'Unnecessary treatments avoided/yr vs SOC\n({start_year}–{end_year} average)',
                  fontsize=11)


# ── Panel C: time series ──────────────────────────────────────────────────────

def plot_unnecessary_ts(dfs, ax, start_year=TS_START, end_year=TS_END):
    """Unnecessary treatments over time for TS_SCENARIOS."""
    for scen in TS_SCENARIOS:
        if scen not in dfs:
            continue
        u     = total_unnecessary(dfs[scen], start_year, end_year)
        label = SCENARIO_LABELS.get(scen, scen).replace('\n', ' ')
        ax.plot(u.index, u, color=TS_COLORS[scen], linewidth=2.5, label=label)

    ax.set_ylabel('Unnecessary treatments per year')
    ax.set_title('(C) Unnecessary treatments over time')
    ax.legend(frameon=False, fontsize=12)
    ax.set_xlim(start_year, end_year)
    ax.set_ylim(bottom=0)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    sc.SIticks(ax, axis='y')


# ── Panel D: summary OT rate bars ────────────────────────────────────────────

def plot_summary_ot_rate(dfs, ax, start_year=BAR_START, end_year=BAR_END):
    """Average OT rate 2027–2040 for TS_SCENARIOS."""
    rates, labels, colors = [], [], []
    for scen in TS_SCENARIOS:
        if scen not in dfs:
            continue
        u = total_unnecessary(dfs[scen], start_year, end_year).mean()
        t = total_treated(dfs[scen],     start_year, end_year).mean()
        rates.append(u / t * 100 if t > 0 else 0)
        labels.append(SCENARIO_LABELS.get(scen, scen))
        colors.append(TS_COLORS[scen])

    x    = np.arange(len(rates))
    bars = ax.bar(x, rates, color=colors, alpha=0.85, edgecolor='white', width=0.65)
    for bar, rate in zip(bars, rates):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.5,
                f'{rate:.0f}%', ha='center', va='bottom', fontsize=12, fontweight='bold')

    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=11)
    ax.set_ylabel('Overtreatment rate (%)')
    ax.set_title(f'(D) Overtreatment rate\n{start_year}–{end_year} average')
    ax.set_ylim(0, 110)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)


# ── Main ──────────────────────────────────────────────────────────────────────

if __name__ == '__main__':

    dfs = load_scenarios(['soc'] + ALL_NON_SOC)

    if 'soc' not in dfs:
        print('SOC results not found. Run run_scenarios.py first.')
        exit()

    set_font(size=18)
    fig = pl.figure(figsize=(24, 16))
    gs  = GridSpec(2, 2, left=0.06, right=0.98, bottom=0.08, top=0.94,
                   wspace=0.30, hspace=0.40)

    plot_pathway_comparison(dfs, fig.add_subplot(gs[0, 0]))
    plot_ranking(dfs,           fig.add_subplot(gs[0, 1]))
    plot_unnecessary_ts(dfs,    fig.add_subplot(gs[1, 0]))
    plot_summary_ot_rate(dfs,   fig.add_subplot(gs[1, 1]))

    pl.savefig(f'{FIGURES_DIR}/fig4_scenario_comparison.png', dpi=200, bbox_inches='tight')
    print(f'Saved {FIGURES_DIR}/fig4_scenario_comparison.png')
