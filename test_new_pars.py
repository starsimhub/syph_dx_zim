"""
Quick test: do old calibration pars sustain syphilis transmission with new model changes?
(2.1: early latent 24mo, 2.2: male chancre visibility 80%, female CS x1.5)
"""
import numpy as np
import stisim as sti
from run_sims import make_sim, load_calib_pars

pars_df = load_calib_pars()
print(f'Loaded {len(pars_df)} parameter sets')

def check_syph_alive(sim):
    if sim is None: return False
    return float(np.sum(sim.results.syph.new_infections[-60:])) > 0

n_test = 10
base = make_sim(verbose=-1)
msim = sti.make_calib_sims(
    calib_pars=pars_df, sim=base, n_parsets=n_test, check_fn=check_syph_alive,
)

n_alive = len(msim.sims)
print(f'\nResult: {n_alive}/{n_test} parameter sets sustain transmission with new model changes')
for i, sim in enumerate(msim.sims):
    prev = sim.results.syph.prevalence[-1]
    ni = np.sum(sim.results.syph.new_infections[-60:])
    print(f'  Sim {i}: prev={prev:.4f} | new_inf(5yr)={ni:.0f}')
