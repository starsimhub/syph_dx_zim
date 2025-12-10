"""
Create HIV model and interventions
"""

# %% Imports and settings
import pandas as pd
import stisim as sti
import starsim as ss
import sciris as sc


def make_diseases(which='all'):

    diseases = sc.autolist()
    if which in ['all', 'hiv']:
        hiv = sti.HIV(
            beta_m2f=0.02,
            eff_condom=0.85,
            init_prev_data=pd.read_csv('data/init_prev_hiv.csv'),
            rel_init_prev=8.,
        )
        diseases.append(hiv)

    if which in ['all', 'stis']:
        syph = sti.Syphilis(
            beta_m2f=0.9,
            beta_m2c=1.,
            eff_condom=0.5,
            rel_trans_latent_half_life=ss.years(.5),
            anc_detection=1.,
            init_prev_data=pd.read_csv('data/init_prev_syph.csv'),
            init_prev_latent_data=pd.read_csv('data/init_prev_latent_syph.csv'),
        )
        gud = sti.GUDPlaceholder(prevalence=0.05)
        diseases.extend([syph, gud])

    connectors = None
    if which == 'all':
        connectors = [sti.hiv_syph(hiv, syph, rel_sus_hiv_syph=2)]

    return diseases, connectors


