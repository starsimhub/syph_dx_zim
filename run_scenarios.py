"""
Run diagnostic scenarios with calibrated parameters.
Saves treatment outcomes for analysis and plotting.

Because syphilis dynamics are marginal with calibrated parameters,
each parameter set is run with multiple seeds and only runs where
syphilis survives are kept (mirroring the calibration's check_fn).
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
from run_sims import make_sim, make_sim_pars

LOCATION = 'zimbabwe'
RESULTS_DIR = 'results'


def check_syph_alive(sim):
    """Check that syphilis didn't die out (same as calibration check_fn)"""
    syph_ni = sim.results.syph.new_infections[-60:]  # Last 5 years
    return float(np.sum(syph_ni)) > 0


def run_scenario(scenario='soc', n_pars=10, seeds_per_par=5, start=1985, stop=2040, do_save=True):
    """
    Run a scenario with multiple calibrated parameter sets.

    For each parameter set, tries multiple seeds and keeps only runs
    where syphilis survives (matching calibration check_fn behavior).
    """

    # Load calibrated parameters, sorted by effective force (strongest first)
    pars_df = sc.loadobj(f'{RESULTS_DIR}/{LOCATION}_pars_all.df')
    pars_df['eff_force'] = pars_df['syph_beta_m2f'] * pars_df['syph_rel_trans_primary'] * (1 - pars_df['syph_eff_condom'])
    pars_df = pars_df.sort_values('eff_force', ascending=False).reset_index(drop=True)
    n_pars = min(n_pars, len(pars_df))
    print(f'Running scenario "{scenario}" with {n_pars} parameter sets x {seeds_per_par} seeds')

    # For each parameter set, create one base sim and init+apply pars once per seed.
    # Deep-copying the UNINITIALIZED base avoids redundant make_sim() calls,
    # while each copy gets its own init() with a unique seed.
    sims = sc.autolist()
    for par_idx in range(n_pars):
        # Create uninitialized base (no sim.init yet)
        base = make_sim(
            dislist='all',
            scenario=scenario,
            seed=1,
            start=start,
            stop=stop,
            verbose=-1,
        )
        calib_pars = pars_df.iloc[par_idx].to_dict()

        for seed in range(1, seeds_per_par + 1):
            sim = sc.dcp(base)
            sim.pars['rand_seed'] = seed
            sim = make_sim_pars(sim, calib_pars)  # Sets pre-init pars, calls init(), sets post-init pars
            sim.par_idx = par_idx
            sim.seed = seed
            sim.scenario = scenario
            sims += sim

        print(f'  Created {seeds_per_par} sims for par_idx {par_idx}/{n_pars} (eff_force={pars_df.iloc[par_idx]["eff_force"]:.3f})')

    print(f'Created {len(sims)} sims, running in parallel...')

    # Run in parallel with progress tracking
    print(f'Running {len(sims)} sims in parallel...')
    sims = ss.parallel(sims, progress_bar=True).sims
    print(f'Completed {len(sims)} simulations')

    # Filter: keep only runs where syphilis survived
    # For each par_idx, keep the first surviving seed
    kept = []
    for par_idx in range(n_pars):
        par_sims = [s for s in sims if s.par_idx == par_idx]
        found = False
        for s in par_sims:
            if check_syph_alive(s):
                kept.append(s)
                found = True
                break
        if not found:
            print(f'  WARNING: par_idx {par_idx} (eff_force={pars_df.iloc[par_idx]["eff_force"]:.3f}) — syphilis died in all {seeds_per_par} seeds, skipping')

    print(f'Kept {len(kept)}/{n_pars} parameter sets with sustained syphilis')

    if do_save and len(kept) > 0:
        save_treatment_outcomes(kept, scenario)

    return kept


def save_treatment_outcomes(sims, scenario):
    """Extract treatment_outcomes analyzer results and save as dataframe"""

    all_rows = []
    for sim in sims:
        par_idx = sim.par_idx
        tx = sim.results.treatment_outcomes
        yearvec = np.array([float(y) for y in sim.t.yearvec])

        for key in tx.keys():
            vals = np.array(tx[key][:], dtype=float)
            for year, val in zip(yearvec, vals):
                all_rows.append(dict(
                    scenario=scenario,
                    par_idx=par_idx,
                    year=year,
                    metric=key,
                    value=val,
                ))

    df = pd.DataFrame(all_rows)
    fname = f'{RESULTS_DIR}/treatment_outcomes_{scenario}.df'
    sc.saveobj(fname, df)
    print(f'Saved treatment outcomes to {fname}')

    return df


if __name__ == '__main__':

    scenarios = ['soc']  # Start with SOC only
    n_pars = 20  # Top 20 parameter sets by effective force
    seeds_per_par = 5  # Try 5 seeds each

    for scenario in scenarios:
        sims = run_scenario(
            scenario=scenario,
            n_pars=n_pars,
            seeds_per_par=seeds_per_par,
            stop=2040,
        )
