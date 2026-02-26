"""
Run diagnostic scenarios with calibrated parameters.
Saves treatment outcomes for analysis and plotting.
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


def run_scenario(scenario='soc', n_pars=10, start=1985, stop=2040, n_agents=10e3, do_save=True):
    """Run a scenario with multiple calibrated parameter sets"""

    # Load calibrated parameters
    pars_df = sc.loadobj(f'{RESULTS_DIR}/{LOCATION}_pars_all.df')
    n_pars = min(n_pars, len(pars_df))
    print(f'Running scenario "{scenario}" with {n_pars} parameter sets')

    sims = sc.autolist()
    for par_idx in range(n_pars):
        sim = make_sim(
            dislist='all',
            scenario=scenario,
            seed=1,  # Use calibration seed to avoid stochastic die-out
            start=start,
            stop=stop,
            verbose=-1,
        )
        # Apply calibrated parameters
        calib_pars = pars_df.iloc[par_idx].to_dict()
        sim = make_sim_pars(sim, calib_pars)
        sim.par_idx = par_idx
        sim.scenario = scenario
        sims += sim

    # Run in parallel
    sims = ss.parallel(sims).sims
    print(f'Completed {len(sims)} simulations')

    if do_save:
        save_treatment_outcomes(sims, scenario)

    return sims


def save_treatment_outcomes(sims, scenario):
    """Extract treatment_outcomes analyzer results and save as dataframe"""

    all_rows = []
    for sim in sims:
        par_idx = sim.par_idx
        tx = sim.results.treatment_outcomes
        yearvec = np.array([float(y) for y in sim.t.yearvec])

        # Get all result keys from the analyzer
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
    n_pars = 10  # Number of calibrated parameter sets to use
    n_agents = 10e3

    for scenario in scenarios:
        sims = run_scenario(
            scenario=scenario,
            n_pars=n_pars,
            n_agents=n_agents,
            stop=2040,
        )
