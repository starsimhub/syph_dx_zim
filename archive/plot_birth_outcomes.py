"""
Plot birth outcomes vs maternal infection duration.

Runs a small set of calibrated sims with the birth_outcome_duration
analyzer, then plots outcome proportions by infection duration bin.
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
import starsim as ss
import pandas as pd
import matplotlib.pyplot as pl
from matplotlib.patches import Patch
import stisim as sti
from run_sims import make_sim
from run_msim import check_syph_alive
from analyzers import birth_outcome_duration
from utils import set_font

LOCATION = 'zimbabwe'
RESULTS_DIR = 'results'


def collect_birth_outcome_data(n_pars=20, seeds_per_par=3, start=1985, stop=2026):
    """Run sims with birth_outcome_duration analyzer and collect events."""
    pars_df = sc.loadobj(f'{RESULTS_DIR}/{LOCATION}_pars_all.df')

    # Sort by effective transmission force (strongest first) for better survival
    pars_df['eff_force'] = pars_df['syph_beta_m2f'] * pars_df['syph_rel_trans_primary'] * (1 - pars_df['syph_eff_condom'])
    pars_df = pars_df.sort_values('eff_force', ascending=False).reset_index(drop=True)
    n_pars = min(n_pars, len(pars_df))
    print(f'Running {n_pars} pars x {seeds_per_par} seeds to collect birth outcome data...')

    sims = sc.autolist()
    for par_idx in range(n_pars):
        calib_pars = pars_df.iloc[par_idx].to_dict()

        base = make_sim(
            scenario='soc',
            seed=1,
            start=start,
            stop=stop,
            verbose=-1,
            analyzers=[birth_outcome_duration()],
        )
        sti.set_sim_pars(base, calib_pars)

        for seed in range(1, seeds_per_par + 1):
            sim = sc.dcp(base)
            sim.pars['rand_seed'] = seed
            sim.par_idx = par_idx
            sims += sim

    sims = ss.parallel(sims).sims

    # Keep first surviving seed per par_idx
    kept = []
    for par_idx in range(n_pars):
        par_sims = [s for s in sims if s.par_idx == par_idx]
        for s in par_sims:
            if check_syph_alive(s):
                kept.append(s)
                break
    print(f'Kept {len(kept)}/{n_pars} parameter sets')

    # Collect all events across sims
    all_events = []
    for sim in kept:
        az = sim.analyzers.get('birth_outcome_duration')
        if az is not None:
            all_events.extend(az.events)

    df = pd.DataFrame(all_events)
    print(f'Collected {len(df)} birth outcome events')
    return df


def plot_birth_outcomes(df, save=True):
    """Plot birth outcome proportions by infection duration."""
    outcome_keys = ['miscarriage', 'nnd', 'stillborn', 'congenital', 'normal']
    outcome_labels = ['Miscarriage', 'Neonatal death', 'Stillbirth', 'Congenital syphilis', 'Normal']
    outcome_colors = ['#984ea3', '#ff7f00', '#555555', '#cb181d', '#4daf4a']

    # Duration bins (years)
    bin_edges = [0, 0.5, 1, 2, 4, 100]
    bin_labels = ['<6m', '6m–1y', '1–2y', '2–4y', '4y+']

    df['dur_bin'] = pd.cut(df['duration'], bins=bin_edges, labels=bin_labels, right=False)

    set_font(size=18)
    fig, axes = pl.subplots(1, 2, figsize=(16, 7), gridspec_kw={'width_ratios': [2, 1]})

    # --- Panel A: Stacked bars of outcome proportions by duration bin ---
    ax = axes[0]
    groups = df.groupby('dur_bin', observed=True)
    x = np.arange(len(bin_labels))
    width = 0.65
    bottom = np.zeros(len(bin_labels))
    n_per_bin = groups.size()

    for oi, (outcome, label, color) in enumerate(zip(outcome_keys, outcome_labels, outcome_colors)):
        if outcome == 'miscarriage':
            continue  # Very rare, skip
        props = []
        for bl in bin_labels:
            if bl in groups.groups:
                g = groups.get_group(bl)
                props.append((g['outcome'] == oi).sum() / len(g) * 100)
            else:
                props.append(0)
        props = np.array(props)
        ax.bar(x, props, width, bottom=bottom, color=color, alpha=0.85,
               edgecolor='white', linewidth=0.5, label=label)

        # Add percentage labels for segments > 5%
        for j in range(len(bin_labels)):
            if props[j] > 5:
                ax.text(x[j], bottom[j] + props[j] / 2, f'{props[j]:.0f}%',
                        ha='center', va='center', fontsize=14, color='white', fontweight='bold')
        bottom += props

    ax.set_xticks(x)
    ax.set_xticklabels(bin_labels)
    ax.set_xlabel('Duration of maternal infection at delivery')
    ax.set_ylabel('Birth outcomes (%)')
    ax.set_title('Birth outcomes by\nmaternal infection duration')
    ax.set_ylim(0, 105)
    ax.legend(frameon=False, fontsize=14, loc='upper right')

    # Add sample sizes
    for j, bl in enumerate(bin_labels):
        n = n_per_bin.get(bl, 0)
        ax.text(x[j], 102, f'n={n:,}', ha='center', va='bottom', fontsize=12, color='grey')

    # --- Panel B: Distribution of infection durations ---
    ax2 = axes[1]
    durations = df['duration'].values
    durations_clipped = np.clip(durations, 0, 15)
    ax2.hist(durations_clipped, bins=30, color='#4a90d9', alpha=0.7, edgecolor='white')
    ax2.set_xlabel('Infection duration at delivery (years)')
    ax2.set_ylabel('Number of births')
    ax2.set_title('Distribution of\ninfection duration')
    ax2.axvline(x=np.median(durations), color='k', linestyle='--', linewidth=1.5,
                label=f'Median: {np.median(durations):.1f}y')
    ax2.legend(frameon=False, fontsize=14)

    pl.tight_layout()
    if save:
        pl.savefig('figures/birth_outcomes_by_duration.png', dpi=200, bbox_inches='tight')
        print('Saved figures/birth_outcomes_by_duration.png')
    pl.show()


if __name__ == '__main__':

    cachefile = f'{RESULTS_DIR}/birth_outcome_events.df'

    # Collect data (or load cached)
    if os.path.exists(cachefile):
        print(f'Loading cached data from {cachefile}')
        df = sc.loadobj(cachefile)
    else:
        df = collect_birth_outcome_data(n_pars=20)
        sc.saveobj(cachefile, df)
        print(f'Cached to {cachefile}')

    if len(df) > 0:
        plot_birth_outcomes(df)
    else:
        print('No birth outcome events collected')
