"""
Run multi-sim with top calibrated parameter sets.

This is the main analysis pipeline — separate from calibration.
Calibration finds the best parameters; this script runs them and
generates all the results needed for figures and tables.

sim.to_df() captures ALL results automatically (744+ columns):
diseases, analyzers (treatment_outcomes, epi_ts, transmission_by_stage,
sw_stats, coinfection_stats), and interventions. No explicit results
list needed — just add a new analyzer and rerun this script.

Usage:
    python run_msim.py              # Run top 200 pars, generate stats
    python run_msim.py --n_pars 50  # Quick run with fewer pars
"""

import os
os.environ.update(
    OMP_NUM_THREADS='1',
    OPENBLAS_NUM_THREADS='1',
    NUMEXPR_NUM_THREADS='1',
    MKL_NUM_THREADS='1',
)

import sciris as sc
import stisim as sti
import pandas as pd
from run_sims import make_sim, load_calib_pars
from utils import percentiles, check_sim_alive

LOCATION = 'zimbabwe'
RESULTS_DIR = 'results'


def _run_one_sim(pars_row, scenario, start, stop):
    """Run a single sim replaying calibrated parameters.

    Uses the rand_seed stored in pars_row, which is the Optuna-suggested seed
    that was applied to the calibration sim via default_build_fn. Falls back
    to seed=1 for legacy pars_df files that predate the reseed fix.
    """
    seed = int(pars_row.get('rand_seed', 1))
    sim = make_sim(scenario=scenario, start=start, stop=stop, verbose=-1, seed=seed)
    sti.set_sim_pars(sim, pars_row)
    sim.init()
    sim.run()
    sim.par_idx = int(pars_row.get('index', 0))
    return sim


def run_msim(n_pars=None, start=1985, stop=2027, scenario='soc', n_workers=None):
    """
    Run all calibrated parameter sets.

    Each parset is replayed with its stored rand_seed (the Optuna-suggested seed
    applied during calibration via default_build_fn). n_pars=None uses all available
    parameter sets.
    """
    pars_df = load_calib_pars(n=n_pars)
    print(f'Running {len(pars_df)} parameter sets (using stored rand_seed per parset)')
    args = [(row.to_dict(), scenario, start, stop) for _, row in pars_df.iterrows()]
    sims = sc.parallelize(_run_one_sim, args, parallelizer='multiprocess', ncpus=n_workers)
    sims = [s for s in sims if check_sim_alive(s)]
    print(f'Kept {len(sims)}/{len(pars_df)} sims (syph+HIV sustained)')
    return sims


def prune_columns(df):
    """
    Drop columns we don't need for any figure or analysis.
    Keeps ~100 columns from ~744. Safe to extend KEEP_PREFIXES if new
    analyzers are added — anything not matching gets dropped.
    """
    KEEP_PREFIXES = [
        # Core epi (calibration fig, Fig 1)
        'time', 'n_alive',
        'hiv.prevalence', 'hiv.n_infected', 'hiv.new_infections', 'hiv.new_deaths',
        'hiv.n_on_art', 'hiv.n_diagnosed',
        'syph.prevalence', 'syph.active_prevalence', 'syph.n_active', 'syph.n_infected',
        'syph.new_infections', 'syph.new_reinfections',
        'syph.new_congenital', 'syph.new_congenital_deaths',
        'syph.new_treated', 'syph.new_treated_unnecessary',
        'syph.detected_pregnant_prevalence', 'syph.delivery_prevalence',
        'syph.new_nnds', 'syph.new_stillborns',

        # Analyzers (Figs 1, 2, S1, S2)
        'epi_ts.',
        'coinfection_stats.',
        'active_coinfection_stats.',
        'transmission_by_stage.',
        'sw_stats.',

        # Treatment outcomes (Figs 2, 3) — keep totals only, drop sex/HIV disaggregation
        'treatment_outcomes.gud_syndromic_treated', 'treatment_outcomes.gud_syndromic_success',
        'treatment_outcomes.gud_syndromic_unnecessary', 'treatment_outcomes.gud_syndromic_failure',
        'treatment_outcomes.anc_screen_treated', 'treatment_outcomes.anc_screen_success',
        'treatment_outcomes.anc_screen_unnecessary', 'treatment_outcomes.anc_screen_failure',
        'treatment_outcomes.kp_screen_treated', 'treatment_outcomes.kp_screen_success',
        'treatment_outcomes.kp_screen_unnecessary', 'treatment_outcomes.kp_screen_failure',
        'treatment_outcomes.plhiv_screen_treated', 'treatment_outcomes.plhiv_screen_success',
        'treatment_outcomes.plhiv_screen_unnecessary', 'treatment_outcomes.plhiv_screen_failure',
        'treatment_outcomes.newborn_treated', 'treatment_outcomes.newborn_success',
        'treatment_outcomes.newborn_unnecessary', 'treatment_outcomes.newborn_failure',
        'treatment_outcomes.secondary_rash_missed',
        'treatment_outcomes.n_active', 'treatment_outcomes.n_infected',
        'treatment_outcomes.gud_syndromic_missed', 'treatment_outcomes.anc_screen_missed',
        'treatment_outcomes.kp_screen_missed', 'treatment_outcomes.plhiv_screen_missed',
        'treatment_outcomes.newborn_missed',
        # Stage at detection
        'treatment_outcomes.gud_syndromic_stage_', 'treatment_outcomes.anc_screen_stage_',
        'treatment_outcomes.kp_screen_stage_', 'treatment_outcomes.plhiv_screen_stage_',
    ]

    keep = [c for c in df.columns if any(c.startswith(p) or c == p for p in KEEP_PREFIXES)]
    dropped = len(df.columns) - len(keep)
    if dropped:
        print(f'  Pruned {dropped} columns → {len(keep)} kept')
    return df[keep]


def _sim_to_df(sim):
    """Extract and prune a single sim's results. Used by save_results for parallel processing."""
    df = sim.to_df(resample='year', use_years=True, sep='.')
    df = prune_columns(df)
    df['par_idx'] = sim.par_idx
    return df


def save_results(sims):
    """
    Generate percentile statistics and save.
    Uses sim.to_df() then prunes to ~100 columns used by plot scripts.
    Parallelises the per-sim to_df/prune step across all available CPUs.
    """
    print('Generating results from sims...')

    dfs = sc.parallelize(_sim_to_df, sims, parallelizer='thread')
    resdf = pd.concat(dfs)
    print(f'  Combined DataFrame: {len(resdf)} rows, {len(resdf.columns)} columns')

    # Generate percentile statistics grouped by year
    cs = resdf.groupby(resdf.timevec).describe(percentiles=percentiles)
    sc.saveobj(f'{RESULTS_DIR}/{LOCATION}_calib_stats_all.df', cs)
    print(f'Saved {RESULTS_DIR}/{LOCATION}_calib_stats_all.df')

    # SW prevalence by age/sex — stored per-sim, not in to_df()
    save_sw_prev(sims)

    return cs


def save_sw_prev(sims):
    """Extract SW prevalence snapshots from sw_prev_snapshot analyzer."""
    all_prev = []
    for sim in sims:
        snap = sim.analyzers.get('sw_prev_snapshot')
        if snap is not None and snap.prev_data:
            for row in snap.prev_data:
                row = dict(row)
                row['par_idx'] = sim.par_idx
                all_prev.append(row)

    if all_prev:
        sw_prev_df = pd.DataFrame(all_prev)
        sc.saveobj(f'{RESULTS_DIR}/sw_prev_df.df', sw_prev_df)
        print(f'Saved {RESULTS_DIR}/sw_prev_df.df ({len(sw_prev_df)} rows)')
    else:
        print('WARNING: No sw_prev data found — sw_prev_df.df not saved')


if __name__ == '__main__':

    n_pars = 200  # 10% of 2000 calibration trials; load_calib_pars caps at available if fewer survive
    stop = 2027

    sims = run_msim(n_pars=n_pars, stop=stop)

    if len(sims) > 0:
        cs = save_results(sims)
    else:
        print('No surviving sims — cannot generate stats')

    print('Done!')
