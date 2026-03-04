"""
Run syphilis-HIV coinfection model
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
FIGURES_DIR = 'figures'


def make_sim(scenario='soc', seed=1, start=1985, stop=2031, verbose=1/12, analyzers=None):

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
        use_migration=False, rand_seed=seed, n_agents=10e3, start=start, stop=stop, verbose=verbose,
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

    seed = 1
    scenario = 'soc'

    sim = make_sim(stop=2031, seed=seed, scenario=scenario)
    sim.run()

    df = sim.to_df(resample='year', use_years=True, sep='.')
    sc.saveobj(f'results/{scenario}_sim.df', df)
