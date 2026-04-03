"""
Run diagnostic scenarios with calibrated parameters.
Saves treatment outcomes for analysis and plotting.

Each scenario replays the same calibrated parameter sets (with their stored
rand_seed) under a different intervention configuration. run_scenario delegates
to run_msim.run_msim() for the execution loop — so scenario='soc' with the same
stop year produces identical sims to run_msim.

Usage:
    python run_scenarios.py              # All 4 scenarios, all pars
    python run_scenarios.py --n_pars 50  # Quick run with fewer pars
"""

import os
os.environ.update(
    OMP_NUM_THREADS='1',
    OPENBLAS_NUM_THREADS='1',
    NUMEXPR_NUM_THREADS='1',
    MKL_NUM_THREADS='1',
)

import sciris as sc
import pandas as pd
from run_msim import run_msim

LOCATION = 'zimbabwe'
RESULTS_DIR = 'results'

# Full 2^4 power set of {gud, anc, kp, plhiv} — 16 scenarios total.
# Names are underscore-joined components; 'soc', 'conf', 'both' are kept as readable aliases.
# Without GUD (8):
#   soc={}  anc={anc}  kp={kp}  plhiv={plhiv}
#   anc_kp  anc_plhiv  kp_plhiv  conf={anc,kp,plhiv}
# With GUD (8):
#   gud  gud_anc  gud_kp  gud_plhiv
#   gud_anc_kp  gud_anc_plhiv  gud_kp_plhiv  both={gud,anc,kp,plhiv}
ALL_SCENARIOS = ['soc', 'gud', 'anc', 'kp', 'plhiv']

SCENARIOS = ALL_SCENARIOS  # default: run everything


def run_scenario(scenario='soc', n_pars=None, start=1985, stop=2041, n_workers=None, do_save=True):
    """
    Run a scenario with top calibrated parameter sets.

    Delegates to run_msim() for the execution loop — same seeds, same filter —
    so differences between scenarios reflect only the intervention.
    """
    sims = run_msim(n_pars=n_pars, start=start, stop=stop, scenario=scenario, n_workers=n_workers)

    if do_save and len(sims) > 0:
        save_treatment_outcomes(sims, scenario)

    return sims


def save_treatment_outcomes(sims, scenario):
    """Extract treatment_outcomes results and save as a wide-format dataframe."""
    dfs = []
    for sim in sims:
        df = sim.results.treatment_outcomes.to_df(resample='year', use_years=True)
        df.index.name = 'year'
        df = df.reset_index()
        if pd.api.types.is_datetime64_any_dtype(df['year']):
            df['year'] = df['year'].dt.year
        df['par_idx'] = sim.par_idx
        df['scenario'] = scenario
        dfs.append(df)

    combined = pd.concat(dfs, ignore_index=True)
    fname = f'{RESULTS_DIR}/treatment_outcomes_{scenario}.df'
    sc.saveobj(fname, combined)
    print(f'Saved treatment outcomes to {fname}')
    return combined


if __name__ == '__main__':

    n_pars = None  # None = all surviving parsets
    stop = 2041

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--scenarios', nargs='+', default=None,
                        help='Subset of scenarios to run (default: all 16)')
    parser.add_argument('--n_pars', type=int, default=None)
    args = parser.parse_args()

    n_pars = args.n_pars
    to_run = args.scenarios if args.scenarios else ALL_SCENARIOS

    n_total = len(to_run)
    for i, scenario in enumerate(to_run):
        print(f'\n[{i+1}/{n_total}] Running scenario: {scenario}')
        run_scenario(scenario=scenario, n_pars=n_pars, stop=stop)
    print(f'\nDone — {n_total}/{n_total} scenarios completed.')
