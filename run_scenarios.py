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
import stisim as sti
import pandas as pd
from run_sims import make_sim
from run_msim import check_syph_alive

LOCATION = 'zimbabwe'
RESULTS_DIR = 'results'


def run_scenario(scenario='soc', n_pars=10, seeds_per_par=5, start=1985, stop=2040, do_save=True):
    """
    Run a scenario with multiple calibrated parameter sets.

    For each parameter set, tries multiple seeds and keeps only runs
    where syphilis survives (matching calibration check_fn behavior).
    """
    pars_df = sc.loadobj(f'{RESULTS_DIR}/{LOCATION}_pars_all.df')
    base = make_sim(scenario=scenario, start=start, stop=stop, verbose=-1)
    msim = sti.make_calib_sims(
        calib_pars=pars_df, sim=base, n_parsets=n_pars,
        seeds_per_par=seeds_per_par, check_fn=check_syph_alive,
    )
    kept = msim.sims
    print(f'Kept {len(kept)} parameter sets with sustained syphilis')

    if do_save and len(kept) > 0:
        save_treatment_outcomes(kept, scenario)

    return kept


def save_treatment_outcomes(sims, scenario):
    """Extract treatment_outcomes analyzer results + congenital outcomes and save as dataframe"""

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

        # Also save congenital syphilis outcomes from the disease results
        syph_res = sim.results.syph
        for key in ['new_congenital', 'new_congenital_deaths']:
            if key in syph_res:
                vals = np.array(syph_res[key][:], dtype=float)
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

    scenarios = ['soc', 'gud', 'conf', 'both', 'cs']
    n_pars = 10  # Top 10 parameter sets by effective force
    seeds_per_par = 1  # Single seed (no variation)

    for scenario in scenarios:
        sims = run_scenario(
            scenario=scenario,
            n_pars=n_pars,
            seeds_per_par=seeds_per_par,
            stop=2041,
        )
