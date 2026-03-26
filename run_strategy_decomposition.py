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

import sciris as sc
import stisim as sti
import pandas as pd
import pylab as pl
from run_sims import make_sim, load_calib_pars
from run_msim import check_syph_alive
from utils import set_font
from run_scenarios import save_treatment_outcomes

LOCATION = 'zimbabwe'
RESULTS_DIR = 'results'


def run_decomposition(n_pars=10, seeds_per_par=3, stop=2041):
    """Run SOC + each strategy in isolation."""
    pars_df = load_calib_pars()

    scenarios = ['soc', 'gud', 'conf_anc', 'conf_fsw', 'conf_plhiv']
    labels = {
        'soc': 'Standard of care',
        'gud': 'GUD POC diagnostic',
        'conf_anc': 'Confirmatory (ANC)',
        'conf_fsw': 'High-risk HIV- dual testing',
        'conf_plhiv': 'Confirmatory (PLHIV)',
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


ADULT_PATHWAYS = ['gud_syndromic', 'anc_screen', 'kp_screen', 'plhiv_screen']


def get_ot_rate_total(df, start_year, end_year):
    """Compute overall adult OT rate as a single percentage."""
    post = df[(df.year >= start_year) & (df.year <= end_year)]
    treated = sum(post[post.metric == f'{pw}_treated'].value.sum() for pw in ADULT_PATHWAYS)
    unnecessary = sum(post[post.metric == f'{pw}_unnecessary'].value.sum() for pw in ADULT_PATHWAYS)
    return unnecessary / max(treated, 1) * 100


def load_and_plot(start_year=2028, end_year=2039):
    """
    Two-panel figure:
      (A) Overall strategy comparison: SOC vs GUD vs confirmatory (all channels)
      (B) Confirmatory test decomposition by delivery channel
    """
    all_scenarios = ['soc', 'gud', 'conf', 'conf_anc', 'conf_fsw', 'conf_plhiv']

    # Load treatment outcomes
    dfs = {}
    for scenario in all_scenarios:
        fname = f'{RESULTS_DIR}/treatment_outcomes_{scenario}.df'
        df = sc.loadobj(fname)
        if len(df) == 0:
            print(f'WARNING: {fname} is empty — no surviving sims for {scenario}. Rerun with more seeds.')
            return None
        dfs[scenario] = df

    # Compute OT rates
    rates = {s: get_ot_rate_total(dfs[s], start_year, end_year) for s in all_scenarios}
    soc_rate = rates['soc']

    # --- Plot ---
    set_font(size=18)
    fig, (ax1, ax2) = pl.subplots(1, 2, figsize=(18, 5), gridspec_kw={'width_ratios': [1, 1.2]})

    # Panel A: Strategy comparison (SOC, GUD, Confirmatory all channels)
    panel_a = [
        ('soc',  'Standard of care',       '#555555'),
        ('gud',  'GUD syphilis detection',  '#e41a1c'),
        ('conf', 'Active syphilis confirmation\nfollowing treponemal screen+', '#377eb8'),
    ]
    for i, (scen, label, color) in enumerate(panel_a):
        ax1.barh(i, rates[scen], color=color, alpha=0.85, edgecolor='white',
                 linewidth=0.5, height=0.6)
        ax1.text(rates[scen] + 0.5, i, f'{rates[scen]:.0f}%',
                 ha='left', va='center', fontsize=16, fontweight='bold')

    ax1.set_yticks(range(len(panel_a)))
    ax1.set_yticklabels([p[1] for p in panel_a])
    ax1.set_xlabel(f'Adult overtreatment rate (%), {start_year}\u2013{end_year}')
    ax1.set_title('(A) Strategy comparison')
    ax1.set_xlim(0, soc_rate * 1.15)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.invert_yaxis()

    # Panel B: Confirmatory test decomposed by delivery channel
    # Show reduction in OT rate vs SOC for each channel in isolation
    combined_reduction = soc_rate - rates['conf']
    channels = [
        ('conf_anc',    'ANC',              '#377eb8'),
        ('conf_fsw',    'HIV- FSW (dual)',  '#ff7f00'),
        ('conf_plhiv',  'PLHIV',            '#4daf4a'),
    ]
    # Sort by reduction (biggest impact first)
    channels.sort(key=lambda x: rates[x[0]])

    # First bar: combined (all channels)
    ax2.barh(0, combined_reduction, color='#377eb8', alpha=0.4, edgecolor='#377eb8',
             linewidth=1.5, height=0.6, linestyle='--')
    ax2.text(combined_reduction + 0.3, 0, f'{combined_reduction:.0f} pp',
             ha='left', va='center', fontsize=16, fontweight='bold', color='#555555')

    # Individual channel bars
    for i, (scen, label, color) in enumerate(channels):
        reduction = soc_rate - rates[scen]
        ax2.barh(i + 1, reduction, color=color, alpha=0.85, edgecolor='white',
                 linewidth=0.5, height=0.6)
        ax2.text(reduction + 0.3, i + 1, f'{reduction:.0f} pp',
                 ha='left', va='center', fontsize=16, fontweight='bold')

    all_labels = ['All channels combined'] + [c[1] for c in channels]
    ax2.set_yticks(range(len(all_labels)))
    ax2.set_yticklabels(all_labels)
    ax2.set_xlabel(f'Reduction in overtreatment rate (percentage points)')
    ax2.set_title('(B) Active syphilis confirmation: impact by delivery channel')
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.invert_yaxis()
    ax2.set_xlim(0, combined_reduction * 1.2)

    pl.tight_layout()
    fig.savefig(f'figures/strategy_decomposition.png', dpi=200, bbox_inches='tight')
    print(f'Saved figures/strategy_decomposition.png')

    # Print summary
    print(f'\nOvertreatment rates ({start_year}\u2013{end_year}):')
    print(f'  SOC:                   {soc_rate:.1f}%')
    print(f'  GUD detection:         {rates["gud"]:.1f}%')
    print(f'  Confirmatory (all):    {rates["conf"]:.1f}%')
    print(f'\nConfirmatory test by channel (reduction vs SOC):')
    for scen, label, _ in channels:
        print(f'  {label:20s}  {soc_rate - rates[scen]:+.1f} pp  →  {rates[scen]:.1f}%')

    return fig


if __name__ == '__main__':

    n_pars = 10
    seeds_per_par = 3

    # Only run the new decomposed scenarios — reuse existing soc, gud results
    to_run = False  # ['conf_anc', 'conf_fsw', 'conf_plhiv']

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
    # if fig is not None:
    #     pl.show()
    print('Done!')
