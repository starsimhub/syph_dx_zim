"""
Run calibrations for the STI model
"""

# Additions to handle numpy multithreading
import os
os.environ.update(
    OMP_NUM_THREADS='1',
    OPENBLAS_NUM_THREADS='1',
    NUMEXPR_NUM_THREADS='1',
    MKL_NUM_THREADS='1',
)

# %% Imports and settings
import numpy as np
import sciris as sc
import stisim as sti
import pandas as pd
from run_sims import make_sim, make_sim_pars


# Constants
LOCATION = 'zimbabwe'
DATA_DIR = 'data'
RESULTS_DIR = 'results'
FIGURES_DIR = 'figures'


# Run settings
TOTAL_TRIALS = 1200
storage = None
do_shrink = True  # Whether to shrink the calibration results


def make_calibration(which='hiv'):

    # Define the calibration parameters
    ckw = dict(suggest_type='suggest_float')  #
    calib_par_dict = dict(
        hiv=dict(
            hiv_beta_m2f=dict(low=0.003, high=0.02, guess=0.008, **ckw),
            hiv_eff_condom=dict(low=0.5, high=0.95, guess=0.75, **ckw),
            hiv_rel_init_prev=dict(low=2, high=8, guess=4, **ckw),
        ),
        syph=dict(
            syph_beta_m2f=dict(low=0.05, high=0.3, guess=0.15, **ckw),
            syph_eff_condom=dict(low=0.2, high=0.7, guess=0.5, **ckw),
            syph_rel_trans_primary=dict(low=3, high=10, guess=5, **ckw),
        ),
        testing=dict(
            rel_symp_test=dict(low=0.5, high=2.0, guess=1.0, **ckw),
            rel_anc_test=dict(low=0.5, high=2.0, guess=1.0, **ckw),
            rel_kp_test=dict(low=0.5, high=2.0, guess=1.0, **ckw),
        ),
    )
    if which != 'all': calib_pars = calib_par_dict[which]
    else:
        calib_pars = sc.objdict()
        for k, v in calib_par_dict.items():
            calib_pars.update(v)

    # Extra results to save
    sres = sc.autolist()
    savedis = ['hiv'] if which == 'hiv' else ['hiv', 'syph']
    for dis in savedis:
        for res in ['prevalence', 'new_infections', 'n_infected']:
            for sk in ['', '_f']:
                sres += dis+'.'+res+sk
    sres += ['hiv.n_diagnosed', 'hiv.n_on_art', 'n_alive']

    if which in ['syph', 'all']:
        sres += [
            'syph.active_prevalence',
            'syph.active_prevalence_f',
            'syph.active_prevalence_m',
            'syph.detected_pregnant_prevalence',
            'syph.n_active',
            'syph.new_treated',
            'syph.new_treated_unnecessary',
            'syph.new_congenital',
            'syph.new_congenital_deaths',
            'coinfection_stats.syph_prev_has_hiv',
            'coinfection_stats.syph_prev_no_hiv',
        ]

    if which == 'syph':
        sres += [
            'hiv.new_deaths',
            'hiv.new_infections',
            'hiv.prevalence_15_49',
        ]

    # Make the sim
    dislist = which if which == 'hiv' else 'all'
    pre_load_calibs = ['hiv'] if which == 'syph' else None
    sim = make_sim(dislist=dislist, pre_load_calibs=pre_load_calibs, verbose=-1, seed=1)
    data = pd.read_csv(f'data/{LOCATION}_{which}_data.csv')

    weights = {
        'syph.n_infected':10.0,
        'syph.new_infections':10.0,
        'syph.active_prevalence':10.0,
    }

    # Reject trials where syphilis dies out
    def check_syph_alive(sim):
        if sim is None:
            return False
        if 'syph' in sim.diseases:
            # Check syphilis is still transmitting in the last 5 years
            syph_ni = sim.results.syph.new_infections[-60:]  # Last 5 years (monthly timesteps)
            if np.sum(syph_ni) == 0:
                return False
        return True

    check_fn = check_syph_alive if which in ['syph', 'all'] else None

    # Make the calibration — use continue_db and keep_db for HPC crash recovery
    calib = sti.Calibration(
        calib_pars=calib_pars,
        extra_results=sres,
        build_fn=make_sim_pars,
        weights=weights,
        sim=sim,
        data=data,
        check_fn=check_fn,
        study_name=f'{LOCATION}_{which}_calibration',
        total_trials=TOTAL_TRIALS,
        die=False, reseed=False, storage=storage, save_results=True,
        continue_db=True, keep_db=True,
    )

    return sim, calib


def run_calibration(calib, which='hiv', do_save=False):

    # Run the calibration
    printstr = f'Running calibration for {which}, {TOTAL_TRIALS} trials'
    sc.heading(printstr)
    calib.calibrate()
    if do_save: sc.saveobj(f'{RESULTS_DIR}/{LOCATION}_calib_{which}.obj', calib)
    print(f'Best pars are {calib.best_pars}')

    return calib


if __name__ == '__main__':

    which = 'all'  # 'hiv', 'syph', or 'all'
    do_run = False
    make_stats = True

    # Run calibration — with continue_db=True, this will resume from any previous run
    if do_run:
        sim, calib = make_calibration(which)
        calib = run_calibration(calib, which)

        print(f'... finished calibration: {which}')
        print(f'Best pars are {calib.best_pars}')

        # Save the results
        # Clean up the database now that results are saved
        calib.remove_db()

        print('Shrinking and saving...')
        if do_shrink:
            sc.saveobj(f'{RESULTS_DIR}/{LOCATION}_calib_{which}_BIG.obj', calib)
            calib = calib.shrink(n_results=int(TOTAL_TRIALS//10))
            sc.saveobj(f'{RESULTS_DIR}/{LOCATION}_calib_{which}.obj', calib)
        else:
            sc.saveobj(f'{RESULTS_DIR}/{LOCATION}_calib_{which}.obj', calib)

        sc.saveobj(f'{RESULTS_DIR}/{LOCATION}_pars_{which}.df', calib.df)

    else:
        calib = sc.loadobj(f'{RESULTS_DIR}/{LOCATION}_calib_{which}.obj')

    if make_stats:
        print('Making stats...')
        from utils import percentiles
        df = calib.resdf
        df_stats = df.groupby(df.time).describe(percentiles=percentiles)
        sc.saveobj(f'{RESULTS_DIR}/{LOCATION}_calib_stats_{which}.df', df_stats)
        par_stats = calib.df.describe(percentiles=[0.05, 0.95])
        sc.saveobj(f'{RESULTS_DIR}/{LOCATION}_par_stats_{which}.df', par_stats)

    print('Done!')


