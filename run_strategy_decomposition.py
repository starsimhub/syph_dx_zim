"""
Decompose the impact of each diagnostic strategy in isolation.

Produces a figure ranking the size of overtreatment reduction for:
  - GUD POC diagnostic (replaces syndromic management)
  - Confirmatory test for ANC positives only
  - Confirmatory test for KP/PLHIV positives only
  - Newborn congenital syphilis testing

Each strategy is compared against SOC baseline.
"""

import os
os.environ.update(
    OMP_NUM_THREADS='1',
    OPENBLAS_NUM_THREADS='1',
    NUMEXPR_NUM_THREADS='1',
    MKL_NUM_THREADS='1',
)

import numpy as np
import sciris as sc
import stisim as sti
import pandas as pd
import pylab as pl
from run_sims import make_sim, load_calib_pars
from run_msim import check_syph_alive
from run_scenarios import save_treatment_outcomes

LOCATION = 'zimbabwe'
RESULTS_DIR = 'results'


def run_decomposition(n_pars=10, seeds_per_par=3, stop=2041):
    """Run SOC + each strategy in isolation."""
    pars_df = load_calib_pars()

    scenarios = ['soc', 'gud', 'conf_anc', 'conf_kp', 'cs']
    labels = {
        'soc': 'Standard of care',
        'gud': 'GUD POC diagnostic',
        'conf_anc': 'Confirmatory (ANC)',
        'conf_kp': 'Confirmatory (KP/PLHIV)',
        'cs': 'Newborn CS testing',
    }

    all_results = {}
    for scenario in scenarios:
        sc.heading(f'Running {labels[scenario]}')
        base = make_sim(scenario=scenario, stop=stop, verbose=-1)
        msim = sti.make_calib_sims(
            calib_pars=pars_df, sim=base, n_parsets=n_pars,
            seeds_per_par=seeds_per_par, check_fn=check_syph_alive,
        )
        print(f'  Kept {len(msim.sims)} sims')
        save_treatment_outcomes(msim.sims, scenario)
        all_results[scenario] = msim.sims

    return all_results, labels


def load_and_plot(scenarios=None, labels=None):
    """Load treatment outcomes and plot strategy comparison."""
    if scenarios is None:
        scenarios = ['soc', 'gud', 'conf_anc', 'conf_kp', 'cs']
    if labels is None:
        labels = {
            'soc': 'Standard of care',
            'gud': 'GUD POC diagnostic',
            'conf_anc': 'Confirmatory (ANC)',
            'conf_kp': 'Confirmatory (KP/PLHIV)',
            'cs': 'Newborn CS testing',
        }

    # Load treatment outcomes
    dfs = {}
    for scenario in scenarios:
        fname = f'{RESULTS_DIR}/treatment_outcomes_{scenario}.df'
        df = sc.loadobj(fname)
        if len(df) == 0:
            print(f'WARNING: {fname} is empty — no surviving sims for {scenario}. Rerun with more seeds.')
            return None
        dfs[scenario] = df

    # Compute overtreatment rate post-intervention (2028-2040)
    post_year = 2028
    ot_rates = {}
    for scenario in scenarios:
        df = dfs[scenario]
        post = df[df.year >= post_year]

        # Sum treated and unnecessary across all pathways
        treated = post[post.metric.str.endswith('_treated')].groupby('par_idx').value.sum()
        unnecessary = post[post.metric.str.endswith('_unnecessary')].groupby('par_idx').value.sum()
        ot_rate = unnecessary / treated.clip(lower=1)
        ot_rates[scenario] = ot_rate

    # Compute reduction in OT rate relative to SOC
    soc_median = ot_rates['soc'].median()
    strategy_scenarios = [s for s in scenarios if s != 'soc']

    strategy_labels = []
    reductions_median = []
    reductions_lo = []
    reductions_hi = []
    for scenario in strategy_scenarios:
        reduction = ot_rates['soc'] - ot_rates[scenario]  # Per par_idx
        reductions_median.append(reduction.median())
        reductions_lo.append(reduction.quantile(0.25))
        reductions_hi.append(reduction.quantile(0.75))
        strategy_labels.append(labels[scenario])

    # Sort by median reduction
    order = np.argsort(reductions_median)[::-1]
    strategy_labels = [strategy_labels[i] for i in order]
    reductions_median = [reductions_median[i] for i in order]
    reductions_lo = [reductions_lo[i] for i in order]
    reductions_hi = [reductions_hi[i] for i in order]

    # Plot
    fig, ax = pl.subplots(figsize=(8, 5))
    y = np.arange(len(strategy_labels))
    xerr = [
        [m - lo for m, lo in zip(reductions_median, reductions_lo)],
        [hi - m for m, hi in zip(reductions_median, reductions_hi)],
    ]
    ax.barh(y, reductions_median, xerr=xerr, color='steelblue', capsize=4, height=0.6)
    ax.set_yticks(y)
    ax.set_yticklabels(strategy_labels)
    ax.set_xlabel('Reduction in overtreatment rate vs SOC')
    ax.set_title(f'Impact of each strategy in isolation\n(SOC overtreatment rate: {soc_median:.0%})')
    ax.axvline(0, color='k', linewidth=0.5)
    ax.invert_yaxis()
    pl.tight_layout()
    fig.savefig(f'figures/strategy_decomposition.png', dpi=150)
    print(f'Saved figures/strategy_decomposition.png')

    return fig


if __name__ == '__main__':

    n_pars = 10
    seeds_per_par = 3

    # Only run the new decomposed scenarios — reuse existing soc, gud, cs results
    to_run = ['conf_anc', 'conf_kp']

    if to_run:
        pars_df = load_calib_pars()
        for scenario in to_run:
            sc.heading(f'Running {scenario}')
            base = make_sim(scenario=scenario, stop=2041, verbose=-1)
            msim = sti.make_calib_sims(
                calib_pars=pars_df, sim=base, n_parsets=n_pars,
                seeds_per_par=seeds_per_par, check_fn=check_syph_alive,
            )
            print(f'  Kept {len(msim.sims)} sims')
            if len(msim.sims) > 0:
                save_treatment_outcomes(msim.sims, scenario)

    fig = load_and_plot()
    if fig is not None:
        pl.show()
    print('Done!')
