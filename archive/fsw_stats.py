"""
FSW sexual behaviour stats from the simulation.
Computes annual partner counts for FSW by inspecting edges at a snapshot year,
plus SW contact rate from the client-side seeking parameters.
"""

import sys
sys.path.insert(0, '/Users/robynstuart/gf/syph_dx_zim')

import numpy as np
import starsim as ss
from run_sims import make_sim

sim = make_sim(scenario='soc', seed=1, start=1985, stop=2026, verbose=0)
sim.run()

net = sim.networks.structuredsexual
ppl = sim.people

print("=" * 60)
print("POPULATION OVERVIEW")
print("=" * 60)
alive = ppl.alive.uids
female_bool = ppl.female[alive]
females = alive[female_bool]
males   = alive[~female_bool]
adult_f = alive[female_bool & (ppl.age[alive] >= 15) & (ppl.age[alive] < 50)]
adult_m = alive[(~female_bool) & (ppl.age[alive] >= 15) & (ppl.age[alive] < 50)]

# FSW and clients (alive)
fsw_alive    = alive[net.fsw[alive]]
client_alive = alive[net.client[alive]]
# Adult (15-49) FSW
fsw_adult    = alive[net.fsw[alive] & female_bool & (ppl.age[alive] >= 15) & (ppl.age[alive] < 50)]
client_adult = alive[net.client[alive] & (~female_bool) & (ppl.age[alive] >= 15) & (ppl.age[alive] < 50)]

print(f"Alive agents:            {len(alive):6d}")
print(f"Adult women 15-49:       {len(adult_f):6d}")
print(f"Adult men 15-49:         {len(adult_m):6d}")
print(f"FSW (alive, any age):    {len(fsw_alive):6d}  ({100*len(fsw_alive)/max(1,len(females)):.1f}% of all females)")
print(f"FSW (15-49):             {len(fsw_adult):6d}  ({100*len(fsw_adult)/max(1,len(adult_f)):.1f}% of adult women)")
print(f"Clients (alive, any):    {len(client_alive):6d}  ({100*len(client_alive)/max(1,len(males)):.1f}% of all males)")
print(f"Clients (15-49):         {len(client_adult):6d}  ({100*len(client_adult)/max(1,len(adult_m)):.1f}% of adult men)")

# Age distribution of FSW
ages_fsw = ppl.age[fsw_adult]
print(f"\n{'='*60}")
print("FSW AGE (15-49 at end of sim)")
print("=" * 60)
print(f"  n={len(ages_fsw)}")
print(f"  mean={np.mean(ages_fsw):.1f}  median={np.median(ages_fsw):.1f}  "
      f"IQR=[{np.percentile(ages_fsw,25):.1f}, {np.percentile(ages_fsw,75):.1f}]  "
      f"min={np.min(ages_fsw):.1f}  max={np.max(ages_fsw):.1f}")

# Concurrent partners (right now in edges)
# SW edges: edge_types 3
print(f"\n{'='*60}")
print("CONCURRENT PARTNERS (from active edges)")
print("=" * 60)
p1 = net.edges.p1
p2 = net.edges.p2
et = net.edges.edge_type

sw_edges_mask = (et == 3)
sw_p1 = p1[sw_edges_mask]  # clients
sw_p2 = p2[sw_edges_mask]  # fsw

# Count SW partners per FSW agent
from collections import Counter
sw_count_per_fsw = Counter(sw_p2.tolist())
# All FSW, including those with zero SW partners right now
sw_partners_now = np.array([sw_count_per_fsw.get(uid, 0) for uid in fsw_adult])
conc_all = net.partners[fsw_adult]

print(f"FSW with ≥1 active edge right now: {np.sum(sw_partners_now > 0)} of {len(fsw_adult)}")
print(f"  SW partners (snapshot):")
print(f"    mean={np.mean(sw_partners_now):.2f}  median={np.median(sw_partners_now):.2f}  "
      f"IQR=[{np.percentile(sw_partners_now,25):.2f}, {np.percentile(sw_partners_now,75):.2f}]  "
      f"max={np.max(sw_partners_now)}")
print(f"  All concurrent partners (all types):")
print(f"    mean={np.mean(conc_all):.2f}  median={np.median(conc_all):.2f}  "
      f"IQR=[{np.percentile(conc_all,25):.2f}, {np.percentile(conc_all,75):.2f}]  "
      f"max={np.max(conc_all):.0f}")

# Annual SW contact rate from parameters
# sw_seeking_rate = probpermonth(1.0) -> mean 1 new client per month
# Each contact is a separate edge (one-time)
# So each FSW's annual SW contacts ≈ n_clients * seeking_rate / n_FSW
dt = sim.dt  # timesteps per year
sw_rate_per_ts = float(net.pars.sw_seeking_rate.to_prob())  # prob per timestep
n_clients_active = len(client_adult)
n_fsw_active = len(fsw_adult)
# Expected client contacts per FSW per year (if uniformly distributed)
contacts_per_fsw_per_year = n_clients_active * sw_rate_per_ts * (1/dt) / max(1, n_fsw_active)

print(f"\n{'='*60}")
print("ANNUAL SW CONTACT RATE (parametric estimate)")
print("=" * 60)
print(f"  sw_seeking_rate (monthly prob): {sw_rate_per_ts:.4f}")
print(f"  Active clients (15-49):         {n_clients_active}")
print(f"  Active FSW (15-49):             {n_fsw_active}")
print(f"  client:FSW ratio:               {n_clients_active/max(1,n_fsw_active):.1f}")
print(f"  Expected SW contacts/FSW/year:  {contacts_per_fsw_per_year:.0f}")

# Lifetime SW partners
lt_sw = net.lifetime_sw_partners[fsw_adult]
ages_at_debut = net.debut[fsw_adult]
career_years  = np.maximum(0, ppl.age[fsw_adult] - ages_at_debut)
annual_rate   = lt_sw / np.maximum(1, career_years)

print(f"\n{'='*60}")
print("LIFETIME SW PARTNERS & CAREER RATE")
print("=" * 60)
print(f"  Lifetime SW partnerships:")
print(f"    mean={np.mean(lt_sw):.1f}  median={np.median(lt_sw):.1f}  "
      f"IQR=[{np.percentile(lt_sw,25):.1f}, {np.percentile(lt_sw,75):.1f}]  max={np.max(lt_sw):.0f}")
print(f"  Annual SW partnership rate (lifetime / career years):")
print(f"    mean={np.mean(annual_rate):.1f}  median={np.median(annual_rate):.1f}  "
      f"IQR=[{np.percentile(annual_rate,25):.1f}, {np.percentile(annual_rate,75):.1f}]")
