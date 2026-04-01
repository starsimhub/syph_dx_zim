"""
Test that run_msim can exactly reproduce calibration results.

Verifies that:
(a) default_build_fn applies rand_seed to each calibration trial (reseed=True fix)
(b) different trials use different seeds (not all seed=1)
(c) replaying a calibrated parset with its stored rand_seed produces an
    identical mismatch to what was recorded during calibration

Usage:
    python test_calibration_reproducibility.py
"""

import numpy as np
import sciris as sc
import stisim as sti

from run_sims import make_sim
from run_calibrations import make_calibration

N_TRIALS   = 10
STUDY_NAME = 'zimbabwe_calibration_test'
DB_NAME    = 'results/zimbabwe_calibration_test.db'


def run_test():
    print('=' * 60)
    print('Test: calibration reproducibility (reseed fix)')
    print('=' * 60)

    # ------------------------------------------------------------------ #
    # Step 1: Run a small calibration                                     #
    # ------------------------------------------------------------------ #
    print(f'\nStep 1: Running calibration ({N_TRIALS} trials, reseed=True)...')

    _, calib = make_calibration()
    calib.run_args.n_trials    = N_TRIALS
    calib.run_args.n_workers   = 1
    calib.run_args.study_name  = STUDY_NAME
    calib.run_args.storage     = f'sqlite:///{DB_NAME}'
    calib.run_args.db_name     = DB_NAME
    calib.run_args.continue_db = False
    calib.run_args.keep_db     = False

    calib.calibrate()

    surviving = calib.df[np.isfinite(calib.df['mismatch'])]
    print(f'  {len(surviving)}/{N_TRIALS} trials produced finite mismatch')

    if len(surviving) == 0:
        print('  SKIP: No surviving parsets — increase N_TRIALS')
        return

    # ------------------------------------------------------------------ #
    # Step 2: Verify seeds are distinct across trials                     #
    # ------------------------------------------------------------------ #
    print('\nStep 2: Checking that trials used different seeds...')
    seeds = calib.df['rand_seed'].values
    print(f'  Seeds: {seeds}')
    n_unique = len(set(seeds))
    assert n_unique > 1, f'All trials used the same seed! Seeds: {seeds}'
    print(f'  {n_unique}/{len(seeds)} unique seeds — PASS')

    # ------------------------------------------------------------------ #
    # Step 3: Replay each surviving parset using its stored rand_seed     #
    # ------------------------------------------------------------------ #
    print('\nStep 3: Replaying surviving parsets with stored rand_seed...')

    n_checked = 0
    for _, row in surviving.iterrows():
        pars_dict   = row.to_dict()
        trial_idx   = int(row['index'])
        stored_seed = int(row['rand_seed'])

        replay_sim = make_sim(scenario='soc', start=1985, stop=2031, verbose=-1, seed=stored_seed)
        sti.set_sim_pars(replay_sim, pars_dict)
        replay_sim.init()
        replay_sim.run()

        replay_mismatch = calib.eval_fn(replay_sim, **calib.eval_kw)
        calib_mismatch  = row['mismatch']

        print(f'  Parset {trial_idx} (seed={stored_seed}): '
              f'calib={calib_mismatch:.6f}  replay={replay_mismatch:.6f}', end='')

        assert np.isclose(replay_mismatch, calib_mismatch, rtol=1e-6), \
            f'MISMATCH DIFFERS: calib={calib_mismatch:.8f}, replay={replay_mismatch:.8f}'

        print('  PASS')
        n_checked += 1

    print(f'\nAll {n_checked} parsets reproduced exactly.')
    print('\nTest PASSED.')


if __name__ == '__main__':
    run_test()
