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
from run_sims import make_sim, _extract_value, _get_disease, _get_network, _get_connector

LOCATION = 'zimbabwe'
RESULTS_DIR = 'results'


def check_syph_alive(sim):
    """Check that syphilis didn't die out (same as calibration check_fn)"""
    syph_ni = sim.results.syph.new_infections[-60:]  # Last 5 years
    return float(np.sum(syph_ni)) > 0


def set_preinit_pars(sim, calib_pars):
    """Set disease/network/connector pars on an uninitialized sim (fast, no init)"""
    for k, pars in calib_pars.items():
        if k in ['rand_seed', 'index', 'mismatch', 'eff_force']:
            continue
        v = _extract_value(pars)
        if v is None:
            continue
        if 'syph_' in k:
            _get_disease(sim, 'syph').pars[k.replace('syph_', '')] = v
        elif 'hiv_' in k:
            _get_disease(sim, 'hiv').pars[k.replace('hiv_', '')] = v
        elif 'nw_' in k:
            nw = _get_network(sim, 'structuredsexual')
            par_name = k.replace('nw_', '')
            if hasattr(nw.pars[par_name], 'set'):
                nw.pars[par_name].set(v)
            else:
                nw.pars[par_name] = v
        elif 'conn_' in k:
            conn = _get_connector(sim, 'hiv_syph')
            par_name = k.replace('conn_', '')
            conn.pars[par_name] = v


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

    # Create UNINITIALIZED sims with pre-init pars set.
    # sim.init() will happen inside ss.parallel workers, utilizing all cores.
    sims = sc.autolist()
    for par_idx in range(n_pars):
        calib_pars = pars_df.iloc[par_idx].to_dict()

        # Extract intervention rel_test values to pass during sim creation
        rel_symp = _extract_value(calib_pars.get('rel_symp_test', 1.0)) or 1.0
        rel_anc = _extract_value(calib_pars.get('rel_anc_test', 1.0)) or 1.0

        # Create one uninitialized base with intervention pars baked in
        base = make_sim(
            dislist='all',
            scenario=scenario,
            seed=1,
            start=start,
            stop=stop,
            verbose=-1,
        )
        # Set pre-init pars (disease, network, connector) — fast, no init
        set_preinit_pars(base, calib_pars)

        # Set intervention rel_test on the uninitialized intervention objects
        for intv in base.pars.interventions:
            if getattr(intv, 'name', '') == 'symp_algo':
                intv.pars['rel_test'] = rel_symp
            elif getattr(intv, 'name', '') == 'anc_screen':
                intv.pars['rel_test'] = rel_anc
            elif getattr(intv, 'name', '') == 'dual_hiv':
                intv.pars['rel_test'] = _extract_value(calib_pars.get('rel_kp_test', 1.0)) or 1.0

        for seed in range(1, seeds_per_par + 1):
            sim = sc.dcp(base)
            sim.pars['rand_seed'] = seed
            sim.par_idx = par_idx
            sim.seed = seed
            sim.scenario = scenario
            sims += sim

    print(f'Created {len(sims)} uninitialized sims, running in parallel (init + run)...')

    # Run in parallel — each worker will call sim.init() then sim.run()
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
