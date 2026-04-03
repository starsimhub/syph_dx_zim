"""
Plot Figure 4: Scenario comparison — impact of novel diagnostics.

Panels:
    A: Per-pathway SOC vs single-use-case stacked bars (correctly treated +
       overtreated), annotated with % OT change
    B: Total unnecessary treatments avoided vs SOC (summed 2027–2040), ranked
    C: Unnecessary treatments over time by pathway — stacked area;
       top sub-panel = SOC, bottom sub-panel = all products ('both')
    D: Cost threshold analysis — break-even POC price for GUD, KP, PLHIV
"""

import sciris as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as pl
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from matplotlib.patches import Patch
from utils import set_font, get_metric

RESULTS_DIR = 'results'
FIGURES_DIR = 'figures'

# ── BIA constants ──────────────────────────────────────────────────────────────
COST_BPG = 3.0   # $ per BPG treatment
POC_MAX  = 5.0   # x-axis ceiling (USD)

BIA_PATHWAYS = [
    ('gud',   'gud_syndromic', 'GUD syndromic',   '#e41a1c'),
    ('kp',    'kp_screen',     'KP dual RDT',     '#984ea3'),
    ('plhiv', 'plhiv_screen',  'PLHIV screening', '#ff7f00'),
]

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

PATHWAY_COLORS = {
    'gud_syndromic': '#e41a1c',
    'anc_screen':    '#377eb8',
    'kp_screen':     '#984ea3',
    'plhiv_screen':  '#ff7f00',
}


# ── Helpers ───────────────────────────────────────────────────────────────────

def fmt_count(val):
    """Format a count as e.g. '1.16M', '441K', '8.1K'."""
    if val >= 1e6:
        return f'{val/1e6:.2f}M'
    elif val >= 1e3:
        return f'{val/1e3:.1f}K'
    return f'{val:.0f}'


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
    """SOC vs single use-case stacked bars for each adult pathway.

    Colors follow PATHWAY_COLORS throughout:
      - solid fill   = correctly treated
      - faded+hatch  = overtreated
    """
    bar_w   = 0.35
    gap     = 0.08   # gap between the two bars in a group
    spacing = 0.25   # extra space between groups

    group_w = 2 * bar_w + gap + spacing
    group_x = np.arange(len(ADULT_PATHWAYS)) * group_w

    max_y = 0  # track for ylim headroom
    for i, pw in enumerate(ADULT_PATHWAYS):
        scen     = PATHWAY_SCENARIOS[pw]
        pw_color = PATHWAY_COLORS[pw]
        x_soc    = group_x[i]
        x_scen   = group_x[i] + bar_w + gap

        for x, key in [(x_soc, 'soc'), (x_scen, scen)]:
            if key not in dfs:
                continue
            s = pathway_metric(dfs[key], pw, 'success',     start_year, end_year)
            u = pathway_metric(dfs[key], pw, 'unnecessary', start_year, end_year)
            # Correctly treated: solid pathway color
            ax.bar(x, s, bar_w, color=pw_color, alpha=0.85)
            # Overtreated: faded pathway color with hatching
            ax.bar(x, u, bar_w, color=pw_color, alpha=0.30, bottom=s,
                   hatch='////', edgecolor=pw_color, linewidth=0.5)
            total = s + u
            if total > 0:
                ax.text(x, total * 1.02, f'{total/1000:.1f}K', ha='center', va='bottom')
            max_y = max(max_y, total)

        # % OT change annotation — positioned at the SOC bar height
        if 'soc' in dfs and scen in dfs:
            s_soc  = pathway_metric(dfs['soc'],  pw, 'success',     start_year, end_year)
            u_soc  = pathway_metric(dfs['soc'],  pw, 'unnecessary', start_year, end_year)
            u_scen = pathway_metric(dfs[scen],   pw, 'unnecessary', start_year, end_year)
            if u_soc > 0:
                pct     = (u_scen - u_soc) / u_soc * 100
                top_soc = s_soc + u_soc
                ax.text(x_scen, top_soc,
                        f'{pct:+.0f}% OT',
                        ha='center', va='bottom',
                        color=pw_color, fontweight='bold')

    tick_x = group_x + bar_w / 2 + gap / 2
    ax.set_xticks(tick_x)
    ax.set_xticklabels([PATHWAY_LABELS[pw] for pw in ADULT_PATHWAYS])
    ax.set_ylabel(f'Annual average ({start_year}–{end_year})')
    ax.set_title('(A) Per-pathway impact\n(left = SOC,  right = use case alone)')
    ax.set_ylim(bottom=0, top=50000)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    sc.SIticks(ax, axis='y')
    # Legend: treatment status only (pathway color is implicit from x-axis grouping)
    ax.legend(handles=[
        Patch(facecolor='#888888', alpha=0.85, label='Correctly treated'),
        Patch(facecolor='#888888', alpha=0.30, hatch='////', edgecolor='#888888',
              label='Overtreated'),
    ], frameon=False, loc='upper left')


# ── Panel B: ranking (two subpanels) ─────────────────────────────────────────

def plot_ranking(dfs, ax, start_year=BAR_START, end_year=BAR_END):
    """Panel B: ranked horizontal bars — total unnecessary treatments avoided vs SOC (2027–2040 sum)."""
    if 'soc' not in dfs:
        ax.text(0.5, 0.5, 'SOC data not available', ha='center', transform=ax.transAxes)
        return

    soc_u = total_unnecessary(dfs['soc'], start_year, end_year).sum()

    rows = []
    for scen in ALL_NON_SOC:
        if scen not in dfs:
            continue
        u       = total_unnecessary(dfs[scen], start_year, end_year).sum()
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
                fmt_count(val), va='center', fontweight='bold')

    ax.set_yticks(y)
    ax.set_yticklabels([r[0] for r in rows])
    ax.invert_yaxis()
    ax.set_title('(B) Unnecessary treatments avoided vs SOC')
    ax.set_xlabel(f'Total unnecessary treatments avoided vs SOC\n({start_year}–{end_year} sum)')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    sc.SIticks(ax, axis='x')


# ── Panel C: stacked area by pathway ─────────────────────────────────────────

def _stacked_pathway_ax(dfs, ax, scen, subtitle, start_year=TS_START, end_year=TS_END):
    """Draw one stacked-area sub-panel for a single scenario."""
    if scen not in dfs:
        ax.text(0.5, 0.5, f'{scen} data not available', ha='center', transform=ax.transAxes)
        return
    years, ys, colors, labels = None, [], [], []
    for pw in ADULT_PATHWAYS:
        s = get_metric(dfs[scen], f'{pw}_unnecessary', start_year, end_year)
        if years is None:
            years = s.index
        ys.append(s.values)
        colors.append(PATHWAY_COLORS[pw])
        labels.append(PATHWAY_LABELS[pw].replace('\n', ' '))

    ax.stackplot(years, ys, colors=colors, labels=labels, alpha=0.85)
    ax.set_title(subtitle)
    ax.set_ylabel('Unnecessary treatments/yr')
    ax.set_xlim(start_year, end_year)
    ax.set_ylim(bottom=0)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    sc.SIticks(ax, axis='y')
    return ax


def plot_unnecessary_stacked(dfs, gs_c, fig, start_year=TS_START, end_year=TS_END):
    """Panel C: stacked area plots — SOC (top sub-panel) and all products (bottom)."""
    inner = GridSpecFromSubplotSpec(2, 1, subplot_spec=gs_c, hspace=0.45)
    ax_top = fig.add_subplot(inner[0])
    ax_bot = fig.add_subplot(inner[1])
    ax_top = _stacked_pathway_ax(dfs, ax_top, 'soc',  '(C) Unnecessary treatments by pathway — SOC',
                                 start_year, end_year)
    _stacked_pathway_ax(dfs, ax_bot, 'both', '(C) Unnecessary treatments by pathway — All products',
                        start_year, end_year)
    # Legend only on top panel
    if ax_top is not None:
        ax_top.legend(frameon=False, loc='upper left')


# ── Panel D: BIA cost threshold ───────────────────────────────────────────────

def _bia_net_savings(poc_price, n_avoided, n_presenters):
    return COST_BPG * n_avoided - poc_price * n_presenters


def _bia_breakeven(n_avoided, n_presenters):
    return COST_BPG * n_avoided / n_presenters if n_presenters > 0 else np.nan


def _load_bia_par_data(scen, col, start=BAR_START, end=BAR_END):
    """Per-parset cumulative counts for one BIA pathway over [start, end]."""
    soc_df  = sc.loadobj(f'{RESULTS_DIR}/treatment_outcomes_soc.df')
    scen_df = sc.loadobj(f'{RESULTS_DIR}/treatment_outcomes_{scen}.df')
    soc_filt  = soc_df [(soc_df.year  >= start) & (soc_df.year  <= end)]
    scen_filt = scen_df[(scen_df.year >= start) & (scen_df.year <= end)]
    rows = []
    for par_idx in sorted(soc_filt.par_idx.unique()):
        s = soc_filt [soc_filt.par_idx  == par_idx]
        g = scen_filt[scen_filt.par_idx == par_idx]
        n_presenters = s[f'{col}_treated'].sum()
        n_treat_poc  = g[f'{col}_treated'].sum()
        rows.append({'par_idx':      par_idx,
                     'n_presenters': n_presenters,
                     'n_treat_poc':  n_treat_poc,
                     'n_avoided':    n_presenters - n_treat_poc})
    return pd.DataFrame(rows)


def plot_bia(ax):
    """Panel D: cost threshold savings curves for GUD, KP, PLHIV pathways."""
    prices = np.linspace(0, POC_MAX, 300)
    ax.axhline(0, color='black', lw=0.8, ls='--', zorder=1)

    for scen, col, label, color in BIA_PATHWAYS:
        par_df = _load_bia_par_data(scen, col)
        sav_mat = np.array([
            _bia_net_savings(prices, row.n_avoided, row.n_presenters)
            for _, row in par_df.iterrows()
        ])
        med = np.median(sav_mat,      axis=0)
        lo  = np.percentile(sav_mat,  5, axis=0)
        hi  = np.percentile(sav_mat, 95, axis=0)

        be_vals = [_bia_breakeven(r.n_avoided, r.n_presenters) for _, r in par_df.iterrows()]
        be_med  = np.median(be_vals)
        be_lo   = np.percentile(be_vals,  5)
        be_hi   = np.percentile(be_vals, 95)

        legend_label = f'{label}  [break-even ${be_med:.2f}, 90% CI ${be_lo:.2f}–${be_hi:.2f}]'
        ax.fill_between(prices, lo / 1e6, hi / 1e6, alpha=0.15, color=color)
        ax.plot(prices, med / 1e6, color=color, lw=2.5, label=legend_label)
        ax.axvline(be_med, color=color, lw=1.0, ls=':', zorder=2)
        ax.scatter([be_med], [0], color=color, s=60, zorder=5, clip_on=False)

    ax.legend(frameon=False, loc='upper right', fontsize='small')
    ax.set_xlabel('POC test price (USD)')
    ax.set_ylabel(f'Cumulative net savings vs SOC\n{BAR_START}–{BAR_END} (USD millions)')
    ax.set_title('(D) Cost threshold analysis: POC diagnostics\nvs syndromic management')
    ax.set_xlim(0, POC_MAX)
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
    plot_ranking(dfs,            fig.add_subplot(gs[0, 1]))
    plot_unnecessary_stacked(dfs, gs[1, 0], fig)
    plot_bia(fig.add_subplot(gs[1, 1]))

    pl.savefig(f'{FIGURES_DIR}/fig4_scenario_comparison.png', dpi=200, bbox_inches='tight')
    print(f'Saved {FIGURES_DIR}/fig4_scenario_comparison.png')
