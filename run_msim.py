"""
Run multi-sim with top calibrated parameter sets.

This is the main analysis pipeline — separate from calibration.
Calibration finds the best parameters; this script runs them and
generates all the results needed for figures and tables.

sim.to_df() captures ALL results automatically (744+ columns):
diseases, analyzers (treatment_outcomes, epi_ts, transmission_by_stage,
sw_stats, coinfection_stats), and interventions. No explicit results
list needed — just add a new analyzer and rerun this script.

Usage:
    python run_msim.py              # Run top 200 pars, generate stats
    python run_msim.py --n_pars 50  # Quick run with fewer pars
"""

import os
os.environ.update(
    OMP_NUM_THREADS='1',
    OPENBLAS_NUM_THREADS='1',
    NUMEXPR_NUM_THREADS='1',
    MKL_NUM_THREADS='1',
)

import argparse
import numpy as np
import sciris as sc
import starsim as ss
import pandas as pd
from run_sims import make_sim, _extract_value
from run_scenarios import set_preinit_pars, check_syph_alive
from utils import percentiles

LOCATION = 'zimbabwe'
RESULTS_DIR = 'results'


def run_msim(n_pars=200, start=1985, stop=2026, scenario='soc'):
    """
    Run top n_pars calibrated parameter sets.
    No seed variation — each par set is a genuinely distinct fit.
    """

    # Load calibrated parameters, sorted by mismatch (best first)
    pars_df = sc.loadobj(f'{RESULTS_DIR}/{LOCATION}_pars_all.df')
    n_pars = min(n_pars, len(pars_df))
    print(f'Running top {n_pars} parameter sets (scenario={scenario}, {start}-{stop})')

    # Create uninitialized sims
    sims = sc.autolist()
    for par_idx in range(n_pars):
        calib_pars = pars_df.iloc[par_idx].to_dict()

        base = make_sim(
            dislist='all',
            scenario=scenario,
            seed=par_idx + 1,
            start=start,
            stop=stop,
            verbose=-1,
        )
        set_preinit_pars(base, calib_pars)

        # Set intervention rel_test pars
        for intv in base.pars.interventions:
            name = getattr(intv, 'name', '')
            if name == 'symp_algo':
                intv.pars['rel_test'] = _extract_value(calib_pars.get('rel_symp_test', 1.0)) or 1.0
            elif name == 'anc_screen':
                intv.pars['rel_test'] = _extract_value(calib_pars.get('rel_anc_test', 1.0)) or 1.0
            elif name == 'dual_hiv':
                intv.pars['rel_test'] = _extract_value(calib_pars.get('rel_kp_test', 1.0)) or 1.0

        base.par_idx = par_idx
        sims += base

    print(f'Created {len(sims)} uninitialized sims, running in parallel...')
    sims = ss.parallel(sims).sims
    print(f'Completed {len(sims)} simulations')

    # Filter: keep only sims where syphilis survived
    kept = [s for s in sims if check_syph_alive(s)]
    n_died = len(sims) - len(kept)
    if n_died:
        print(f'  Dropped {n_died}/{len(sims)} sims where syphilis died out')
    print(f'Kept {len(kept)} sims')

    return kept


def save_results(sims):
    """
    Generate percentile statistics and save.
    Uses sim.to_df() which captures ALL results automatically (~744 columns).
    Only saves the stats (percentile bands), not the raw per-sim DataFrame.
    """
    print('Generating results from sims...')

    dfs = sc.autolist()
    for i, sim in enumerate(sims):
        df = sim.to_df(resample='year', use_years=True, sep='.')
        df['par_idx'] = sim.par_idx
        dfs += df
        if (i + 1) % 50 == 0:
            print(f'  Processed {i + 1}/{len(sims)} sims')

    resdf = pd.concat(dfs)
    print(f'  Combined DataFrame: {len(resdf)} rows, {len(resdf.columns)} columns')

    # Generate percentile statistics grouped by year
    cs = resdf.groupby(resdf.time).describe(percentiles=percentiles)
    sc.saveobj(f'{RESULTS_DIR}/{LOCATION}_calib_stats_all.df', cs)
    print(f'Saved {RESULTS_DIR}/{LOCATION}_calib_stats_all.df')

    # SW prevalence by age/sex — stored per-sim, not in to_df()
    save_sw_prev(sims)

    return cs


def save_sw_prev(sims):
    """Extract SW prevalence snapshots — these are stored per-sim, not in to_df()"""
    all_prev = []
    for sim in sims:
        sw = sim.analyzers.get('sw_stats')
        if sw is not None and hasattr(sw, 'prev_data') and sw.prev_data is not None:
            for row in sw.prev_data:
                row = dict(row)
                row['par_idx'] = sim.par_idx
                all_prev.append(row)

    if all_prev:
        sw_prev_df = pd.DataFrame(all_prev)
        sc.saveobj(f'{RESULTS_DIR}/sw_prev_df.df', sw_prev_df)
        print(f'Saved {RESULTS_DIR}/sw_prev_df.df ({len(sw_prev_df)} rows)')


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--n_pars', type=int, default=200)
    parser.add_argument('--stop', type=int, default=2026)
    args = parser.parse_args()

    sims = run_msim(n_pars=args.n_pars, stop=args.stop)

    if len(sims) > 0:
        cs = save_results(sims)
    else:
        print('No surviving sims — cannot generate stats')

    print('Done!')
