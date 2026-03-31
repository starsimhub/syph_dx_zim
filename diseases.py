"""
Create HIV model and interventions
"""

# %% Imports and settings
import pandas as pd
import stisim as sti
import starsim as ss


def make_diseases():

    hiv = sti.HIV(
        beta_m2f=0.02,
        eff_condom=0.85,
        init_prev_data=pd.read_csv('data/init_prev_hiv.csv'),
        rel_init_prev=8.,
    )

    syph = sti.Syphilis(
        beta_m2f=0.15,
        beta_m2c=0.075,
        eff_condom=0.5,
        rel_trans_primary=5,                      # High: primary drives ~50-60% of transmission
        rel_trans_secondary=1,                    # Moderate: ~25-30% of transmission
        rel_trans_latent=0.1,                     # Low baseline, decays exponentially
        rel_trans_latent_half_life=ss.months(6),  # Faster decay than default 1yr
        dur_early=ss.uniform(ss.months(22), ss.months(26)),  # WHO definition: early latent = first 24 months (was 12-14)
        p_symp_primary=[0.3, 0.8],               # Chancre visibility: 30% female, 80% male (was 50% male; reflects heterosexual external anatomy)
        anc_detection=1.,
        rel_init_prev=0.2,
        init_prev_data=pd.read_csv('data/init_prev_syph.csv'),
        # init_prev_latent_data=pd.read_csv('data/init_prev_latent_syph.csv'),
    )
    syph.store_sw = True
    gud = sti.GUDPlaceholder(prevalence=0.01)

    diseases = [hiv, syph, gud]
    connectors = [sti.hiv_syph(hiv, syph)]  # Defaults: rel_sus_hiv_syph=2.67, rel_sus_syph_hiv=1; overridden by calibration

    return diseases, connectors
