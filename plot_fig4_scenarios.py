"""
Plot Figure 4: Scenario comparison — impact of novel diagnostics.

Layout: 2 × 3

    A (0,0): Per-pathway SOC vs single-use-case stacked bars (correctly treated +
             overtreated), annotated with % OT change
    B (0,1): Unnecessary treatments over time — stacked area;
             top sub-panel = SOC, bottom sub-panel = all products ('both')
    C (0,2): Total unnecessary treatments avoided vs SOC (summed 2027–2040), ranked
    D (1,0): Unnecessary treatments avoided per test performed, by use case
    E (1,1): Net cost per 1,000 tests vs cost ratio (test cost / treatment cost),
             one line per use case
    F (1,2): Heatmap grid — net savings per test (x: test cost, y: treatment cost),
             2 × 2 sub-panels, one per use case
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
    ('gud',   'gud_syndromic', 'GUD POC test',             '#e41a1c'),
    ('anc',   'anc_screen',    'POC NT diagnostic — ANC',  '#377eb8'),
    ('kp',    'kp_screen',     'POC NT diagnostic — KP',   '#984ea3'),
    ('plhiv', 'plhiv_screen',  'POC NT diagnostic — PLHIV','#ff7f00'),
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

        # Pre-compute % OT change for the scenario bar label
        pct_str = ''
        if 'soc' in dfs and scen in dfs:
            u_soc = pathway_metric(dfs['soc'], pw, 'unnecessary', start_year, end_year)
            u_sc  = pathway_metric(dfs[scen],  pw, 'unnecessary', start_year, end_year)
            if u_soc > 0:
                pct_str = f'\n{(u_sc - u_soc) / u_soc * 100:+.0f}% OT'

        for x, key in [(x_soc, 'soc'), (x_scen, scen)]:
            if key not in dfs:
                continue
            s = pathway_metric(dfs[key], pw, 'success',     start_year, end_year)
            u = pathway_metric(dfs[key], pw, 'unnecessary', start_year, end_year)
            ax.bar(x, s, bar_w, color=pw_color, alpha=0.85)
            ax.bar(x, u, bar_w, color=pw_color, alpha=0.30, bottom=s,
                   hatch='////', edgecolor=pw_color, linewidth=0.5)
            total = s + u
            if total > 0:
                suffix = pct_str if key == scen else ''
                ax.text(x, total * 1.02, f'{total/1000:.1f}K{suffix}',
                        ha='center', va='bottom', fontsize=11,
                        color=pw_color if suffix else 'black', fontweight='bold' if suffix else 'normal')
            max_y = max(max_y, total)

    # Per-bar x-tick labels: GUD\nSOC, GUD\nPOC, ANC\nSOC, etc.
    short = {'gud_syndromic': 'GUD', 'anc_screen': 'ANC',
             'kp_screen': 'KP', 'plhiv_screen': 'PLHIV'}
    tick_xs, tick_lbs = [], []
    for i, pw in enumerate(ADULT_PATHWAYS):
        tick_xs.append(group_x[i])
        tick_xs.append(group_x[i] + bar_w + gap)
        tick_lbs.append(f'{short[pw]}\nSOC')
        tick_lbs.append(f'{short[pw]}\nPOC')
    ax.set_xticks(tick_xs)
    ax.set_xticklabels(tick_lbs, fontsize=14)
    ax.tick_params(axis='x', length=0)

    ax.set_ylabel(f'Annual average ({start_year}–{end_year})')
    ax.set_title(f'(C) Average annual unnecessary treatments\n{BAR_START}–{BAR_END}')
    ax.set_ylim(bottom=0, top=62000)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    sc.SIticks(ax, axis='y')
    # Legend: treatment status only (pathway color is implicit from x-axis grouping)
    ax.legend(handles=[
        Patch(facecolor='#888888', alpha=0.85, label='Correctly treated'),
        Patch(facecolor='#888888', alpha=0.30, hatch='////', edgecolor='#888888',
              label='Overtreated'),
    ], frameon=False, loc='upper left', fontsize=12)


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
    ax.set_title(f'(D) Cumulative unnecessary treatments\navoided vs SOC, {BAR_START}–{BAR_END}')
    ax.set_xlabel(f'Total unnecessary treatments avoided vs SOC\n({start_year}–{end_year} sum)')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    sc.SIticks(ax, axis='x')


# ── Panel C: stacked area by pathway ─────────────────────────────────────────

def _stacked_pathway_ax(dfs, ax, scen, sublabel, start_year=TS_START, end_year=TS_END,
                         hide_xaxis=False):
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
    # Sub-label centred at top of sub-panel
    ax.text(0.5, 0.95, sublabel, transform=ax.transAxes,
            ha='center', va='top', fontsize=15, fontweight='bold', color='#444444')
    ax.set_ylabel('')
    ax.set_xlim(start_year, end_year)
    ax.set_ylim(bottom=0)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    sc.SIticks(ax, axis='y')
    if hide_xaxis:
        ax.tick_params(axis='x', labelbottom=False)
        ax.spines['bottom'].set_visible(False)
    else:
        ax.xaxis.set_major_formatter(pl.FuncFormatter(lambda x, _: f'{int(x)}'))
    return ax


def plot_unnecessary_stacked(dfs, gs_c, fig, start_year=TS_START, end_year=TS_END):
    """Panel A: stacked area plots — SOC (top) and all products (bottom), unified panel."""
    inner = GridSpecFromSubplotSpec(2, 1, subplot_spec=gs_c, hspace=0)
    ax_top = fig.add_subplot(inner[0])
    ax_bot = fig.add_subplot(inner[1], sharex=ax_top, sharey=ax_top)
    ax_top = _stacked_pathway_ax(dfs, ax_top, 'soc',  'SOC', start_year, end_year,
                                  hide_xaxis=True)
    _stacked_pathway_ax(dfs, ax_bot, 'both', 'All products', start_year, end_year)

    if ax_top is None:
        return

    # Unified look: remove touching spines
    ax_top.spines['bottom'].set_visible(False)
    ax_bot.spines['top'].set_visible(False)

    # Single panel title — now in row 0, col 0
    ax_top.set_title(f'(A) Unnecessary treatments\n{BAR_START}–{BAR_END}')

    # Single y-label centred across both axes via fig.text
    pos_top = ax_top.get_position()
    pos_bot = ax_bot.get_position()
    mid_y   = (pos_top.y1 + pos_bot.y0) / 2
    ax_top.set_ylabel('')
    ax_bot.set_ylabel('')
    fig.text(pos_bot.x0 - 0.045, mid_y, 'Unnecessary treatments/yr',
             ha='center', va='center', rotation='vertical', fontsize=14)

    # Legend inside — upper right of top panel
    handles, labels = ax_top.get_legend_handles_labels()
    ax_top.legend(handles, labels, frameon=False, loc='upper left',
                  fontsize=12, ncol=1)


# ── Panel D: unnecessary treatments avoided per test ──────────────────────────

def _load_bia_par_data(scen, col, start=BAR_START, end=BAR_END):
    """Per-parset cumulative counts for one use case over [start, end].

    Returns a DataFrame with columns: par_idx, n_presenters, n_treat_poc, n_avoided,
    avoided_per_test (= n_avoided / n_presenters).
    """
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
        n_avoided    = n_presenters - n_treat_poc
        rows.append({'par_idx':          par_idx,
                     'n_presenters':     n_presenters,
                     'n_treat_poc':      n_treat_poc,
                     'n_avoided':        n_avoided,
                     'avoided_per_test': n_avoided / n_presenters if n_presenters > 0 else np.nan})
    return pd.DataFrame(rows)


def plot_bia(ax):
    """Panel D: unnecessary treatments avoided per test performed, by use case."""
    labels, medians, lo_errs, hi_errs, colors = [], [], [], [], []

    for scen, col, label, color in BIA_PATHWAYS:
        par_df = _load_bia_par_data(scen, col)
        vals   = par_df['avoided_per_test'].dropna().values
        med    = np.median(vals)
        lo     = np.percentile(vals,  5)
        hi     = np.percentile(vals, 95)
        labels.append(label)
        medians.append(med)
        lo_errs.append(med - lo)
        hi_errs.append(hi - med)
        colors.append(color)

    x = np.arange(len(labels))
    ax.bar(x, medians, color=colors, alpha=0.85, width=0.6)
    ax.errorbar(x, medians,
                yerr=[lo_errs, hi_errs],
                fmt='none', color='#333333', capsize=5, lw=1.5, zorder=5)

    # Value labels above each bar
    for xi, med, hi_err in zip(x, medians, hi_errs):
        ax.text(xi, med + hi_err + 0.01, f'{med:.2f}',
                ha='center', va='bottom', fontsize=13, fontweight='bold')

    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel('Unnecessary treatments avoided\nper test performed')
    ax.set_title(f'(E) Unnecessary treatments avoided\nper test performed ({BAR_START}–{BAR_END})')
    ax.set_ylim(0, 1.05)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)


# ── Panel E: cost curve ───────────────────────────────────────────────────────

def plot_cost_curve(ax):
    """Panel E: net cost per 1,000 tests vs cost ratio (test cost / treatment cost).

    Net cost = 1000 × C_treat × (r − avoided_per_test), where r = test_cost / C_treat.
    Break-even at r = avoided_per_test. Lines below zero are cost-saving.
    """
    r = np.linspace(0, 1.5, 400)
    ax.axhline(0, color='black', lw=1.0, ls='--', zorder=3)

    for scen, col, label, color in BIA_PATHWAYS:
        par_df = _load_bia_par_data(scen, col)
        vals   = par_df['avoided_per_test'].dropna().values
        if len(vals) == 0:
            continue
        med = np.median(vals)
        lo  = np.percentile(vals,  5)
        hi  = np.percentile(vals, 95)

        net_med = 1000 * COST_BPG * (r - med)
        net_lo  = 1000 * COST_BPG * (r - hi)   # hi avoided_per_test → lower (better) net cost
        net_hi  = 1000 * COST_BPG * (r - lo)

        label_clean = label.replace('\n', ' ')
        ax.plot(r, net_med, color=color, lw=2.5, label=f'{label_clean}  (r* = {med:.2f})')
        ax.fill_between(r, net_lo, net_hi, color=color, alpha=0.15)
        ax.axvline(med, color=color, lw=0.8, ls=':', zorder=2)

    ax.set_xlabel('Cost ratio  (test cost / treatment cost)')
    ax.set_ylabel('Net cost per 1,000 tests (USD)')
    ax.set_title('(F) Net cost vs cost ratio\n(break-even where line crosses zero)')
    ax.set_xlim(0, 1.5)
    ax.legend(frameon=False, fontsize=12, loc='upper left',
              bbox_to_anchor=(0, 1), borderaxespad=0)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)


# ── Panel F: heatmap grid ─────────────────────────────────────────────────────

def plot_heatmap_grid(gs_f, fig, panel_label=True):
    """Panel F: 2×2 grid of heatmaps — net savings per test performed.

    x-axis: test cost ($0–$5)
    y-axis: treatment cost ($0–$10)
    fill:   net savings per test = avoided_per_test × treatment_cost − test_cost
            (positive = cost-saving, negative = net cost)
    """
    from matplotlib.colors import TwoSlopeNorm

    test_costs  = np.linspace(0, 5,  100)
    treat_costs = np.linspace(0, 10, 100)
    T, C = np.meshgrid(test_costs, treat_costs)

    inner = GridSpecFromSubplotSpec(2, 2, subplot_spec=gs_f, hspace=0.35, wspace=0.30)
    axes  = [fig.add_subplot(inner[r, c]) for r in range(2) for c in range(2)]

    vmin, vmax = -5.0, 5.0   # net savings per test, $
    norm = TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)
    cmap = 'RdYlGn'

    for ax, (scen, col, label, _) in zip(axes, BIA_PATHWAYS):
        par_df = _load_bia_par_data(scen, col)
        vals   = par_df['avoided_per_test'].dropna().values
        if len(vals) == 0:
            ax.text(0.5, 0.5, 'No data', ha='center', transform=ax.transAxes)
            continue
        med = np.median(vals)

        net = med * C - T          # net savings per test = avoided_per_test × C_treat − C_test
        net = np.clip(net, vmin, vmax)

        im = ax.pcolormesh(T, C, net, cmap=cmap, norm=norm, shading='auto')
        ax.contour(T, C, net, levels=[0], colors='black', linewidths=1.2, linestyles='--')
        ax.set_title(label.replace('\n', ' '), color='black', fontsize=13)
        ax.set_xlabel('Test cost ($)')
        ax.set_ylabel('Treatment cost ($)')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    # Shared colorbar on the right of the last axis
    fig.colorbar(im, ax=axes, label='Net savings per test ($)', shrink=0.6, pad=0.02)

    if panel_label:
        fig.text((axes[0].get_position().x0 + axes[1].get_position().x1) / 2,
                 axes[0].get_position().y1 + 0.02,
                 '(F) Net savings per test ($)',
                 ha='center', va='bottom', fontsize=18, fontweight='bold')


# ── Panel F: OT rate + correctly treated over time ────────────────────────────

def plot_ot_and_correct(dfs, gs_f, fig, start_year=TS_START, end_year=TS_END):
    """Panel F: two sub-panels — (top) overtreatment rate over time by scenario;
    (bottom) correctly treated over time by scenario."""
    inner   = GridSpecFromSubplotSpec(2, 1, subplot_spec=gs_f, hspace=0)
    ax_top  = fig.add_subplot(inner[0])
    ax_bot  = fig.add_subplot(inner[1], sharex=ax_top)

    scen_order = ['soc', 'gud', 'anc', 'kp', 'plhiv', 'both']
    scen_labels = {k: v.replace('\n', ' ') for k, v in SCENARIO_LABELS.items()}

    for scen in scen_order:
        if scen not in dfs:
            continue
        color = TS_COLORS.get(scen, '#888888')
        lw    = 2.5 if scen in ('soc', 'both') else 1.5
        ls    = '--' if scen == 'soc' else '-'

        # Total unnecessary and success across all adult pathways
        unnec = total_unnecessary(dfs[scen], start_year, end_year)
        succ  = None
        for pw in ADULT_PATHWAYS:
            s = get_metric(dfs[scen], f'{pw}_success', start_year, end_year)
            succ = s if succ is None else succ.add(s, fill_value=0)

        total = unnec.add(succ, fill_value=0)
        ot_rate = (unnec / total * 100).where(total > 0)

        ax_top.plot(unnec.index, ot_rate, color=color, lw=lw, ls=ls,
                    label=scen_labels[scen])
        ax_bot.plot(succ.index, succ, color=color, lw=lw, ls=ls)

    ax_top.set_ylabel('OT rate (%)')
    ax_top.set_ylim(0, 105)
    ax_top.set_title('(B) Overtreatment rate & correctly treated\nby scenario')
    ax_top.tick_params(axis='x', labelbottom=False)
    ax_top.spines['top'].set_visible(False)
    ax_top.spines['right'].set_visible(False)
    ax_top.spines['bottom'].set_visible(False)
    handles, labels = ax_top.get_legend_handles_labels()
    ax_bot.legend(handles, labels, frameon=False, fontsize=12, loc='upper center',
                  bbox_to_anchor=(0.5, -0.25), ncol=3)

    ax_bot.set_ylabel('Correctly treated/yr')
    ax_bot.set_xlim(start_year, end_year)
    ax_bot.set_ylim(bottom=0)
    ax_bot.xaxis.set_major_formatter(pl.FuncFormatter(lambda x, _: f'{int(x)}'))
    ax_bot.spines['top'].set_visible(False)
    ax_bot.spines['right'].set_visible(False)
    sc.SIticks(ax_bot, axis='y')


# ── Main ──────────────────────────────────────────────────────────────────────

if __name__ == '__main__':

    dfs = load_scenarios(['soc'] + ALL_NON_SOC)

    if 'soc' not in dfs:
        print('SOC results not found. Run run_scenarios.py first.')
        exit()

    set_font(size=16)

    # ── Figure 4: 2×3 layout ──
    fig4 = pl.figure(figsize=(22, 13))
    gs4  = GridSpec(2, 3, left=0.07, right=0.97, bottom=0.10, top=0.94,
                    wspace=0.38, hspace=0.50)

    plot_unnecessary_stacked(dfs, gs4[0, 0], fig4)              # A — (0,0)
    plot_ot_and_correct(dfs,      gs4[0, 1], fig4)              # B — (0,1)
    plot_pathway_comparison(dfs,  fig4.add_subplot(gs4[0, 2]))  # C — (0,2)
    plot_ranking(dfs,             fig4.add_subplot(gs4[1, 0]))  # D — (1,0)
    plot_bia(fig4.add_subplot(gs4[1, 1]))                       # E — (1,1)
    plot_cost_curve(fig4.add_subplot(gs4[1, 2]))                # F — (1,2)

    fig4.savefig(f'{FIGURES_DIR}/fig4_scenario_comparison.png', dpi=200, bbox_inches='tight')
    print(f'Saved {FIGURES_DIR}/fig4_scenario_comparison.png')
    pl.close(fig4)

    # ── Figure 5: heatmap grid ──
    fig5 = pl.figure(figsize=(12, 6))
    gs5  = GridSpec(1, 1, left=0.08, right=0.88, bottom=0.08, top=0.96)
    plot_heatmap_grid(gs5[0, 0], fig5, panel_label=False)

    fig5.savefig(f'{FIGURES_DIR}/fig5_heatmaps.png', dpi=200, bbox_inches='tight')
    print(f'Saved {FIGURES_DIR}/fig5_heatmaps.png')
