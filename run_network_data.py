"""
Run sims and extract network structure data for the supplementary network figure.
Network structure is largely input-driven, so only a few sims are needed.
"""

# %% Imports and settings
import sciris as sc
import numpy as np
import starsim as ss
import stisim as sti
import stisim as sti
from run_sims import make_sim
from analyzers import NetworkSnapshot

LOCATION = 'zimbabwe'
RESULTS_DIR = 'results'


def run_network_sims(n_sims=3, par_idx=0):
    """Run a small number of sims with network analyzers and calibrated params"""

    # Load calibration parameters
    pars_df = sc.loadobj(f'{RESULTS_DIR}/{LOCATION}_pars_all.df')
    calib_pars = pars_df.iloc[par_idx].to_dict()
    print(f'Loaded calibrated parameters (index {par_idx})')

    sims = sc.autolist()
    for seed in range(1, n_sims + 1):
        network_analyzers = [
            sti.partner_age_diff(year=2020),
            NetworkSnapshot(year=2020),
        ]
        sim = make_sim(
            seed=seed,
            stop=2026,
            analyzers=network_analyzers,
        )
        sti.set_sim_pars(sim, calib_pars)
        sim.init()
        sim.par_idx = par_idx
        sims += sim

    sims = ss.parallel(sims).sims
    return sims


def extract_network_data(sims):
    """Extract network data from completed sims and save"""
    data = sc.objdict()

    # Panel A: Lifetime partner distributions (debuted agents only, aggregate across sims)
    lp_f = []
    lp_m = []
    for sim in sims:
        snap = sim.analyzers['network_snapshot']
        lp_f.extend(list(snap.lifetime_partners_data['Female']))
        lp_m.extend(list(snap.lifetime_partners_data['Male']))
    data.lifetime_partners_f = np.array(lp_f)
    data.lifetime_partners_m = np.array(lp_m)

    # Panel B: Age differences (aggregate across sims)
    age_diffs = {}
    for sim in sims:
        pad = sim.analyzers['partner_age_diff']
        for key, vals in pad.age_diffs.items():
            age_diffs.setdefault(key, []).extend(list(vals))
    data.age_diffs = {k: np.array(v) for k, v in age_diffs.items()}

    # Panels C, D, E: From NetworkSnapshot (use first sim as representative)
    snap = sims[0].analyzers['network_snapshot']
    data.risk_group_data = snap.risk_group_data
    data.debut_data = snap.debut_data
    data.rel_dur_data = snap.rel_dur_data

    # Panel E: Average partnership-by-age across all sims
    all_stable = []
    all_casual = []
    for sim in sims:
        s = sim.analyzers['network_snapshot']
        all_stable.append(s.partnership_by_age['prop_stable'])
        all_casual.append(s.partnership_by_age['prop_casual'])
    data.partnership_by_age = dict(
        age_bins=snap.partnership_by_age['age_bins'],
        prop_stable=np.nanmean(all_stable, axis=0),
        prop_casual=np.nanmean(all_casual, axis=0),
    )

    sc.saveobj(f'{RESULTS_DIR}/network_data.obj', data)
    print(f'Saved network data to {RESULTS_DIR}/network_data.obj')
    return data


if __name__ == '__main__':
    sims = run_network_sims(n_sims=3, par_idx=0)
    extract_network_data(sims)
