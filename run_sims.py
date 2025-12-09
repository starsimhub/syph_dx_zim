"""
Run syphilis-HIV coinfection model
"""

# %% Imports and settings
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


def make_sim_pars(sim, calib_pars):
    if not sim.initialized: sim.init()
    hiv = sim.diseases.hiv
    nw = sim.networks.structuredsexual
    if 'syph' in sim.diseases: syph = sim.diseases.syph

    # Apply the calibration parameters
    for k, pars in calib_pars.items():  # Loop over the calibration parameters
        if k == 'rand_seed':
            sim.pars.rand_seed = v
            continue

        elif k in ['index', 'mismatch']:
            continue

        if isinstance(pars, dict):
            v = pars['value']
        elif sc.isnumber(pars):
            v = pars
        else:
            raise NotImplementedError(f'Parameter {k} not recognized')

        if 'syph_' in k:  # Syphilis parameters
            k = k.replace('syph_', '')  # Strip off indentifying part of parameter name
            syph.pars[k] = v
        elif 'hiv_' in k:  # HIV parameters
            k = k.replace('hiv_', '')  # Strip off indentifying part of parameter name
            hiv.pars[k] = v
        elif 'nw_' in k:  # Network parameters
            k = k.replace('nw_', '')  # As above
            if 'pair_form' in k:
                nw.pars[k].set(v)
            else:
                nw.pars[k] = v
        else:
            raise NotImplementedError(f'Parameter {k} not recognized')

    return sim


def make_sim(dislist='all', scenario='soc', seed=1, start=1985, stop=2031, verbose=1/12, analyzers=None, pre_load_calibs=None, par_idx=0):

    # Network
    sexual = sti.StructuredSexual(
        prop_f0=0.6,
        prop_f2=0.05,
        prop_m0=0.5,
        f1_conc=0.05,
        f2_conc=0.25,
        m1_conc=0.15,
        m2_conc=0.3,
        p_pair_form=0.6,
        condom_data=pd.read_csv(f'data/condom_use.csv'),
    )
    maternal = ss.MaternalNet()
    networks = [sexual, maternal]

    # Diseases
    diseases, connectors = make_diseases(which=dislist)

    # Interventions
    interventions = make_interventions(which=dislist, scenario=scenario)

    # Analyzers
    analyzers = make_analyzers(which=dislist, extra_analyzers=analyzers)

    simpars = dict(
        use_migration=True, rand_seed=seed, rel_death=0.99,
        n_agents=10e3, start=start, stop=stop, verbose=verbose,
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
    if pre_load_calibs is not None:
        calib_folder = 'results'
        calib_pars = {}
        for disease in pre_load_calibs:
            pars_df = sc.loadobj(f'{calib_folder}/{LOCATION}_pars_{disease}.df')
            print(f'Loaded calibration parameters for {disease} from {calib_folder}/{LOCATION}_pars_{disease}.df')
            best_pars = pars_df.iloc[par_idx].to_dict()
            print(best_pars)
            calib_pars.update(best_pars)
        sim.init()
        sim = make_sim_pars(sim, calib_pars)
        print(f'Using calibration parameters for scenario {scenario} and index {par_idx}')

        # Display pars from sim
        for disease in pre_load_calibs:
            print(f'\n{disease.upper()} parameters:')
            for k, v in sim.diseases[disease].pars.items():
                if k in best_pars:
                    print(f'  {k}: {v} (calibrated)')
    return sim


if __name__ == '__main__':

    # Set up and run

    # SETTINGS
    debug = False
    seed = 1
    do_save = True
    do_run = True
    do_plot = True
    use_calib = False
    scenario = 'soc'

    to_run = [
        'run_hiv',
        # 'run_all',
        # 'run_msim',
    ]

    if 'run_hiv' in to_run or 'run_stis' in to_run or 'run_all' in to_run:
        dislist = 'all' if 'run_all' in to_run else 'hiv' if 'run_hiv' in to_run else 'stis'
        if do_run:
            sim = make_sim(dislist=dislist, stop=2040, seed=seed, scenario=scenario, pre_load_calibs=None)
            print(f'Running sim for diseases: {dislist}')
            print('Initializing sim...')
            if not sim.initialized: sim.init()
            print('Diseases in sim:', sim.diseases.keys())
            print('Analyzers in sim:', sim.analyzers.keys())
            sim.run()

            if do_save:
                sc.saveobj(f'results/{scenario}_sim_{dislist}.obj', sim)
                df = sim.to_df(resample='year', use_years=True, sep='.')
                sc.saveobj(f'results/{scenario}_sim_{dislist}.df', df)
        else:
            df = sc.loadobj(f'results/{scenario}_sim_{dislist}.df')

        if do_plot:
            from plot_sims import plot_sims
            df.index = df['timevec']
            plot_sims(df, dislist=dislist, which='single', start_year=1985, title=f'{dislist}_plots')

