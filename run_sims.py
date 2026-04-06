"""
Core sim-building functions for the syphilis diagnostics Zimbabwe analysis.

Provides:
  - make_sim()        — construct a fully configured Sim for a given scenario
  - load_calib_pars() — load calibrated parameter sets from disk
  - LEGACY_PAR_RENAME — column rename map for old underscore-style calib cols
"""

# %% Imports and settings
import numpy as np
import sciris as sc
import pandas as pd
import starsim as ss
import stisim as sti

# From this repo
from diseases import make_diseases
from interventions import make_interventions
from analyzers import make_analyzers

# Constants
LOCATION = 'zimbabwe'
DATA_DIR = 'data'
RESULTS_DIR = 'results'

# Rename map: old underscore-format calibration columns → new dot notation
LEGACY_PAR_RENAME = {
    'rel_symp_test': 'symp_algo.rel_test',
    'rel_anc_test':  'anc_screen.rel_test',
    'rel_kp_test':   'dual_hiv.rel_test',
}


def load_calib_pars(path=None, sort_by_force=True, n=None):
    """
    Load calibrated parameters, renaming legacy column names to dot notation.

    If sort_by_force=True, sorts by effective syphilis transmission force
    (strongest first) so the top-N parameter sets are most likely to sustain
    syphilis dynamics.
    """
    if path is None:
        path = f'{RESULTS_DIR}/{LOCATION}_pars_all.df'
    df = sc.loadobj(path)
    df = df.rename(columns=LEGACY_PAR_RENAME)
    if sort_by_force:
        beta_col = next((c for c in df.columns if c in ['syph.beta_m2f', 'syph_beta_m2f']), None)
        rtp_col = next((c for c in df.columns if c in ['syph.rel_trans_primary', 'syph_rel_trans_primary']), None)
        eff_col = next((c for c in df.columns if c in ['syph.eff_condom', 'syph_eff_condom']), None)
        if all(c is not None for c in [beta_col, rtp_col, eff_col]):
            df['eff_force'] = df[beta_col] * df[rtp_col] * (1 - df[eff_col])
            df = df.sort_values('eff_force', ascending=False).reset_index(drop=True)
    if n is not None:
        df = df.head(n)
    return df


def make_sim(scenario='soc', seed=1, start=1985, stop=2027, verbose=1/12, analyzers=None):

    # Network
    sexual = sti.StructuredSexual(
        prop_f0=0.6,
        prop_f2=0.1,  # 60% LR, 30% MR, 10% HR
        prop_m0=0.5,
        prop_m2=0.1,  # 50% LR, 50% MR, 10% HR
        f1_conc=0.15,
        f2_conc=0.25,
        m1_conc=0.15,
        m2_conc=0.5,
        p_pair_form=0.5,
        condom_data=pd.read_csv(f'data/condom_use.csv'),
        fsw_shares=ss.bernoulli(p=0.10),
        client_shares=ss.bernoulli(p=0.20),
    )
    maternal = ss.MaternalNet()
    networks = [sexual, maternal]

    # Diseases
    diseases, connectors = make_diseases()

    # Interventions
    interventions = make_interventions(scenario=scenario, rel_symp_test=1.2, rel_anc_test=1.1)

    # Analyzers
    analyzers = make_analyzers(extra_analyzers=analyzers)

    simpars = dict(
        use_migration=False, rand_seed=seed, n_agents=10e3, age_scale=1000, start=start, stop=stop, verbose=verbose,
    )

    sim = sti.Sim(
        pars=simpars,
        datafolder='data/',
        demographics=LOCATION.lower(),
        diseases=diseases,
        networks=networks,
        connectors=connectors,
        interventions=interventions,
        analyzers=analyzers,
    )
    sim.scenario = scenario

    print('Created sim')
    return sim


if __name__ == '__main__':
    # Quick smoke-test: build and run a single SOC sim to 2027
    sim = make_sim(stop=2027, seed=1, scenario='soc')
    sim.run()
    print(f'Final syphilis prevalence: {sim.results.syph.active_prevalence[-1]:.4f}')
