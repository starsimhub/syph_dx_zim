"""
Cost-impact curve: net cost per 1,000 tests vs cost ratio (test cost / treatment cost).

For each use case, the net cost per 1,000 tests is:

    net_cost = 1000 × C_treat × (r − avoided_per_test)

where r = test_cost / C_treat is the cost ratio and avoided_per_test is the
fraction of tests that prevent an unnecessary treatment. The break-even cost
ratio is r* = avoided_per_test, i.e. cost-saving whenever the test costs less
than avoided_per_test × C_treat.

Usage:
    python plot_bia_costcurve.py
"""

import sciris as sc
import numpy as np
import matplotlib.pyplot as pl
from utils import set_font

RESULTS_DIR = 'results'
FIGURES_DIR = 'figures'

COST_BPG = 3.0   # $ per BPG treatment course

BAR_START = 2027
BAR_END   = 2040

SCALE = 1000  # "per N tests"

USE_CASES = [
    ('gud',   'gud_syndromic', 'GUD syndromic',   '#e41a1c'),
    ('anc',   'anc_screen',    'ANC screening',   '#377eb8'),
    ('kp',    'kp_screen',     'KP dual RDT',     '#984ea3'),
    ('plhiv', 'plhiv_screen',  'PLHIV screening', '#ff7f00'),
]


def load_avoided_per_test(scen, col, start=BAR_START, end=BAR_END):
    """Per-parset avoided_per_test = n_avoided / n_presenters over [start, end]."""
    soc_df  = sc.loadobj(f'{RESULTS_DIR}/treatment_outcomes_soc.df')
    scen_df = sc.loadobj(f'{RESULTS_DIR}/treatment_outcomes_{scen}.df')
    soc_f  = soc_df [(soc_df.year  >= start) & (soc_df.year  <= end)]
    scen_f = scen_df[(scen_df.year >= start) & (scen_df.year <= end)]
    vals = []
    for par_idx in sorted(soc_f.par_idx.unique()):
        s = soc_f [soc_f.par_idx  == par_idx]
        g = scen_f[scen_f.par_idx == par_idx]
        n_pres   = s[f'{col}_treated'].sum()
        n_poc    = g[f'{col}_treated'].sum()
        if n_pres > 0:
            vals.append((n_pres - n_poc) / n_pres)
    return np.array(vals)


if __name__ == '__main__':
    set_font(size=16)
    fig, ax = pl.subplots(figsize=(8, 6))

    r = np.linspace(0, 1.5, 400)  # cost ratio axis

    for scen, col, label, color in USE_CASES:
        vals = load_avoided_per_test(scen, col)
        if len(vals) == 0:
            continue
        med = np.median(vals)
        lo  = np.percentile(vals,  5)
        hi  = np.percentile(vals, 95)

        net_med = SCALE * COST_BPG * (r - med)
        net_lo  = SCALE * COST_BPG * (r - hi)   # hi avoided → lower (better) net cost
        net_hi  = SCALE * COST_BPG * (r - lo)

        ax.plot(r, net_med, color=color, lw=2.5, label=f'{label}  (break-even r = {med:.2f})')
        ax.fill_between(r, net_lo, net_hi, color=color, alpha=0.15)
        ax.axvline(med, color=color, lw=0.8, ls=':', zorder=2)

    ax.axhline(0, color='black', lw=1.0, ls='--', zorder=3)
    ax.set_xlabel('Cost ratio  (test cost / treatment cost)')
    ax.set_ylabel(f'Net cost impact per {SCALE:,} tests (USD)')
    ax.set_title('Cost-saving if test cost < avoided treatments × treatment cost')
    ax.set_xlim(0, 1.5)
    ax.legend(frameon=False, fontsize=13)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    pl.tight_layout()
    outpath = f'{FIGURES_DIR}/bia_costcurve.png'
    pl.savefig(outpath, dpi=200, bbox_inches='tight')
    print(f'Saved {outpath}')
