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
debug = False  # If True, this will do smaller runs that can be run locally for debugging
n_trials = [100, 2][debug]  # How many trials to run for calibration
n_workers = [100, 1][debug]    # How many cores to use
storage = None
do_shrink = True  # Whether to shrink the calibration results
study_name = 'syph_calibration'


def make_calibration(which='hiv', n_trials=None, n_workers=None):

    # Define the calibration parameters
    ckw = dict(suggest_type='suggest_float')  #
    calib_par_dict = dict(
        hiv=dict(
            hiv_beta_m2f=dict(low=0.008, high=0.02, guess=0.012, **ckw),
            hiv_eff_condom=dict(low=0.5, high=0.95, guess=0.75, **ckw),
            nw_prop_f0 = dict(low=0.55, high=0.9, guess=0.85, **ckw),
            nw_prop_m0 = dict(low=0.50, high=0.9, guess=0.81, **ckw),
            nw_f1_conc = dict(low=0.01, high=0.2, guess=0.01, **ckw),
            nw_m1_conc = dict(low=0.01, high=0.2, guess=0.01, **ckw),
            nw_p_pair_form = dict(low=0.4, high=0.9, guess=0.5, **ckw),
        ),
        syph=dict(
            syph_beta_m2f=dict(low=0.02, high=0.25, guess=0.08, **ckw),
        ),
    )
    calib_pars = calib_par_dict[which]

    # Extra results to save
    sres = sc.autolist()
    savedis = ['hiv'] if which == 'hiv' else ['hiv', 'syph']
    for dis in savedis:
        for res in ['prevalence', 'new_infections', 'n_infected']:
            for sk in ['', '_f']:
                sres += dis+'.'+res+sk
    sres += ['hiv.n_diagnosed', 'hiv.n_on_art', 'n_alive']

    # Make the sim
    dislist = which if which == 'hiv' else 'all'
    pre_load_calibs = ['hiv'] if which == 'syph' else None
    sim = make_sim(dislist=dislist, pre_load_calibs=pre_load_calibs, verbose=-1, seed=1)
    data = pd.read_csv(f'data/{LOCATION}_{which}_data.csv')

    weights = dict(
    )

    # Make the calibration
    calib = sti.Calibration(
        calib_pars=calib_pars,
        extra_results=sres,
        build_fn=make_sim_pars,
        weights=weights,
        sim=sim,
        data=data,
        study_name=f'{LOCATION}_{which}_calibration',
        total_trials=n_trials, n_workers=n_workers,
        die=False, reseed=False, storage=storage, save_results=True,
    )

    return sim, calib


def run_calibration(calib, which='hiv', n_trials=None, do_save=False):

    # Run the calibration
    printstr = f'Running calibration for {which}, {n_trials} trials'
    sc.heading(printstr)
    calib.calibrate()
    if do_save: sc.saveobj(f'{RESULTS_DIR}/{LOCATION}_calib_{which}.obj', calib)
    print(f'Best pars are {calib.best_pars}')

    return calib


if __name__ == '__main__':

    load_partial = False
    which = 'hiv'  # 'hiv' or 'all'
    do_run = True
    do_plot = False
    make_stats = True  # Whether to make stats

    # Run calibration
    if do_run:
        sc.heading(f'Running calibration: {which}')
        sim, calib = make_calibration(which, n_trials=n_trials, n_workers=n_workers)

        if load_partial:
            # Load a partially-run calibration study
            import optuna as op
            print(calib.run_args.study_name)
            study = op.load_study(storage=calib.run_args.storage, study_name=calib.run_args.study_name)
            # calib.run_args.continue_db = True
            # calib.calibrate()
            output = study.optimize(calib.run_trial, n_trials=2)
            calib.best_pars = sc.objdict(study.best_params)
            calib.parse_study(study)
            print('Best pars:', calib.best_pars)

            # Tidy up
            calib.calibrated = True
            if not calib.run_args.keep_db:
                calib.remove_db()

        else:
            calib = run_calibration(calib, which, n_trials=n_trials)

        print(f'... finished calibration: {which}')
        print(f'Best pars are {calib.best_pars}')

        # Save the results
        print('Shrinking and saving...')
        if do_shrink:
            sc.saveobj(f'{RESULTS_DIR}/{LOCATION}_calib_{which}_BIG.obj', calib)
            calib = calib.shrink(n_results=int(n_trials//10))  # Save 5% best results
            sc.saveobj(f'{RESULTS_DIR}/{LOCATION}_calib_{which}.obj', calib)
        else:
            sc.saveobj(f'{RESULTS_DIR}/{LOCATION}_calib_{which}.obj', calib)

        # Save the parameter dataframe
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

    if do_plot:
        from plot_sims import plot_calibrations
        plot_calibrations(dislist=which, start_year=1990)

    print('Done!')


