"""
Budget impact analysis — GUD, KP, and PLHIV POC diagnostics.

Break-even threshold: maximum POC test price at which each scenario
is cost-saving vs SOC, given:
  - BPG (benzathine penicillin G): $3 per treatment
  - Consultation cost: $2 per presentation (both scenarios, cancels out)
  - Drug + test costs only (no CEA)

Cost model per pathway:
  SOC:  n_presenters × $3 BPG
  POC:  n_presenters × cost_poc + n_treated_poc × $3 BPG

Net savings = $3 × n_avoided − cost_poc × n_presenters
  (consultation $2 identical in both arms → cancels)

Break-even: cost_poc* = $3 × n_avoided / n_presenters
"""

import sciris as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as pl
from utils import set_font

RESULTS_DIR  = 'results'
FIGURES_DIR  = 'figures'

COST_BPG = 3.0   # $ per BPG treatment
POC_MAX  = 5.0   # x-axis ceiling

BAR_START = 2027
BAR_END   = 2040

# Pathway config: (scenario file, treatment_outcomes column, display label, color)
PATHWAYS = [
    ('gud',   'gud_syndromic', 'GUD syndromic',  '#e41a1c'),
    ('kp',    'kp_screen',     'KP dual RDT',    '#984ea3'),
    ('plhiv', 'plhiv_screen',  'PLHIV screening','#ff7f00'),
]


# ── Data ─────────────────────────────────────────────────────────────────────

def load_par_data(scen, col, start=BAR_START, end=BAR_END):
    """Per-parset cumulative stats for one pathway over [start, end]."""
    soc_df  = sc.loadobj(f'{RESULTS_DIR}/treatment_outcomes_soc.df')
    scen_df = sc.loadobj(f'{RESULTS_DIR}/treatment_outcomes_{scen}.df')

    soc  = soc_df [(soc_df.year  >= start) & (soc_df.year  <= end)]
    scen = scen_df[(scen_df.year >= start) & (scen_df.year <= end)]

    rows = []
    for par_idx in sorted(soc.par_idx.unique()):
        s = soc [soc.par_idx  == par_idx]
        g = scen[scen.par_idx == par_idx]

        n_presenters = s[f'{col}_treated'].sum()
        n_treat_poc  = g[f'{col}_treated'].sum()
        n_avoided    = n_presenters - n_treat_poc

        rows.append({
            'par_idx':      par_idx,
            'n_presenters': n_presenters,
            'n_treat_poc':  n_treat_poc,
            'n_avoided':    n_avoided,
        })

    return pd.DataFrame(rows)


# ── Cost helpers ─────────────────────────────────────────────────────────────

def net_savings(poc_price, n_avoided, n_presenters):
    return COST_BPG * n_avoided - poc_price * n_presenters


def breakeven_price(n_avoided, n_presenters):
    return COST_BPG * n_avoided / n_presenters if n_presenters > 0 else np.nan


# ── Plot ─────────────────────────────────────────────────────────────────────

def plot_savings_curves(pathway_data, ax):
    """One savings curve per pathway, all on the same axes."""
    prices = np.linspace(0, POC_MAX, 300)

    ax.axhline(0, color='black', lw=0.8, ls='--', zorder=1)

    # Track break-even x positions to space annotations vertically
    for label, color, par_df in pathway_data:
        sav_mat = np.array([
            net_savings(prices, row.n_avoided, row.n_presenters)
            for _, row in par_df.iterrows()
        ])
        med = np.median(sav_mat,     axis=0)
        lo  = np.percentile(sav_mat,  5, axis=0)
        hi  = np.percentile(sav_mat, 95, axis=0)

        be_vals = [breakeven_price(r.n_avoided, r.n_presenters) for _, r in par_df.iterrows()]
        be_med  = np.median(be_vals)
        be_lo   = np.percentile(be_vals,  5)
        be_hi   = np.percentile(be_vals, 95)

        legend_label = f'{label}  [break-even ${be_med:.2f}, 90% CI ${be_lo:.2f}–${be_hi:.2f}]'
        ax.fill_between(prices, lo / 1e6, hi / 1e6, alpha=0.15, color=color)
        ax.plot(prices, med / 1e6, color=color, lw=2.5, label=legend_label)
        ax.axvline(be_med, color=color, lw=1.0, ls=':', zorder=2)
        ax.scatter([be_med], [0], color=color, s=60, zorder=5, clip_on=False)

    ax.legend(frameon=False, loc='upper right')
    ax.set_xlabel('POC test price (USD)')
    ax.set_ylabel(f'Cumulative net savings vs SOC\n{BAR_START}–{BAR_END} (USD millions)')
    ax.set_title('Cost threshold analysis: POC diagnostics\ncompared to syndromic management')
    ax.set_xlim(0, POC_MAX)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)


# ── Main ─────────────────────────────────────────────────────────────────────

if __name__ == '__main__':

    pathway_data = []
    for scen, col, label, color in PATHWAYS:
        par_df = load_par_data(scen, col)
        be = [breakeven_price(r.n_avoided, r.n_presenters) for _, r in par_df.iterrows()]
        print(f'{label:20s}  presenters={np.median(par_df.n_presenters):>8,.0f}'
              f'  avoided={np.median(par_df.n_avoided):>8,.0f}'
              f'  break-even=${np.median(be):.2f}'
              f'  (90% CI ${np.percentile(be,5):.2f}–${np.percentile(be,95):.2f})')
        pathway_data.append((label, color, par_df))

    set_font(size=18)
    fig, ax = pl.subplots(1, 1, figsize=(10, 7))
    fig.subplots_adjust(left=0.12, right=0.97, bottom=0.12, top=0.90)

    plot_savings_curves(pathway_data, ax)

    pl.savefig(f'{FIGURES_DIR}/fig_bia.png', dpi=200, bbox_inches='tight')
    print(f'\nSaved {FIGURES_DIR}/fig_bia.png')
