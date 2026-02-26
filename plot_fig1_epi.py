"""
Plot syphilis and HIV epidemiology in Zimbabwe: Fig 1 in manuscript.

Panels:
    A: Active syphilis prevalence by HIV status + ZIMPHIA
    B: Syphilis transmission by stage (3 bars)
    C: Active syphilis prevalence by age and sex + ZIMPHIA
    D: Syphilis infections by sex work status (%)
    E: HIV infections by sex work status (%)
    F: HIV prevalence by age and sex + ZIMPHIA
"""

import sciris as sc
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as pl
from matplotlib.gridspec import GridSpec
from utils import set_font


# ZIMPHIA 2015-2016 data
ZIMPHIA_SYPH = dict(
    active_hiv_pos=0.029, active_hiv_neg=0.004,
    active_f=0.010, active_m=0.006,
)
# ZIMPHIA 2015-2016 HIV prevalence by age/sex (from published ZIMPHIA report)
ZIMPHIA_HIV = {
    # age_group: (female, male)
    '15-20': (0.040, 0.025),
    '20-25': (0.077, 0.034),
    '25-30': (0.137, 0.077),
    '30-35': (0.207, 0.148),
    '35-50': (0.233, 0.207),
    '50-65': (0.130, 0.142),
}

# Colors — sex-disaggregated
F_COLOR = '#d46e9c'
M_COLOR = '#4a90d9'
F_COLOR_LIGHT = '#f0a3c4'
M_COLOR_LIGHT = '#a3c4e8'

# Colors — HIV-disaggregated
HIV_POS_COLOR = '#e41a1c'
HIV_NEG_COLOR = '#377eb8'


def plot_active_syph_by_hiv(cs, ax, start_year=1995):
    """Panel A: Active syphilis prevalence by HIV status"""
    years = cs.index[cs.index >= start_year]

    for col, label, color in [
        ('active_coinfection_stats.syph_prev_has_hiv', 'HIV+', HIV_POS_COLOR),
        ('active_coinfection_stats.syph_prev_no_hiv', 'HIV\u2212', HIV_NEG_COLOR),
    ]:
        med = cs.loc[years, (col, '50%')] * 100
        lo = cs.loc[years, (col, '10%')] * 100
        hi = cs.loc[years, (col, '90%')] * 100
        ax.fill_between(years, lo, hi, alpha=0.15, color=color)
        ax.plot(years, med, color=color, label=label, linewidth=2)

    # ZIMPHIA data points at 2016
    ax.scatter([2016], [ZIMPHIA_SYPH['active_hiv_pos'] * 100], color=HIV_POS_COLOR,
               marker='D', s=80, zorder=5, edgecolors='k', linewidths=0.5, label='ZIMPHIA 2016')
    ax.scatter([2016], [ZIMPHIA_SYPH['active_hiv_neg'] * 100], color=HIV_NEG_COLOR,
               marker='D', s=80, zorder=5, edgecolors='k', linewidths=0.5)

    ax.legend(frameon=False, fontsize=12, loc='upper left')
    ax.set_ylabel('Prevalence (%)')
    ax.set_title('Active syphilis prevalence\nby HIV status')
    ax.set_ylim(bottom=0)
    ax.set_xlim(start_year, 2025)


def plot_transmission_by_stage_bars(cs, ax):
    """Panel B: 3 bars showing % of sexual transmission by stage"""
    stages = ['primary', 'secondary', 'early']
    labels = ['Primary', 'Secondary', 'Early\nlatent']
    colors = ['#e41a1c', '#ff7f00', '#984ea3']

    # Get cumulative totals across all years
    total = 0
    vals = {}
    for stage in stages + ['late', 'tertiary']:
        col = f'transmission_by_stage.new_sex_{stage}'
        if col in cs.columns.get_level_values(0):
            v = cs[(col, '50%')].sum()
            vals[stage] = v
            total += v

    pcts = [vals.get(s, 0) / total * 100 for s in stages]

    bars = ax.bar(labels, pcts, color=colors, alpha=0.85, edgecolor='white', linewidth=0.5)
    for bar, pct in zip(bars, pcts):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 1,
                f'{pct:.0f}%', ha='center', va='bottom', fontsize=14, fontweight='bold')

    ax.set_ylabel('Share of sexual\ntransmissions (%)')
    ax.set_title('Syphilis transmission\nby disease stage')
    ax.set_ylim(0, 80)
    ax.set_xlabel('')


def plot_active_prev_by_age(epi_df, ax):
    """Panel C: Active syphilis prevalence by age and sex + ZIMPHIA"""
    thisdf = epi_df.loc[(epi_df.disease == 'syph') &
                        (epi_df.age != '0-15') &
                        (epi_df.age != '65+')].copy()
    thisdf['active_prevalence'] *= 100
    sns.barplot(data=thisdf, x='age', y='active_prevalence', hue='sex',
                ax=ax, palette=[F_COLOR, M_COLOR], alpha=0.85)

    # ZIMPHIA reference lines
    ax.axhline(y=ZIMPHIA_SYPH['active_f'] * 100, color=F_COLOR, linestyle='--',
               linewidth=1.5, alpha=0.6, label='ZIMPHIA F')
    ax.axhline(y=ZIMPHIA_SYPH['active_m'] * 100, color=M_COLOR, linestyle='--',
               linewidth=1.5, alpha=0.6, label='ZIMPHIA M')

    ax.set_title('Active syphilis prevalence\nby age')
    ax.set_ylabel('Prevalence (%)')
    ax.set_xlabel('Age group')
    ax.set_ylim(bottom=0)
    ax.legend(frameon=False, fontsize=11)


def plot_infections_by_sw_pct(sw_df, disease, ax, start_year=2000, end_year=2019):
    """Panel D/E: Infections by sex work status as percentages"""
    groups = {'fsw': 'FSW', 'client': 'Client', 'non_fsw': 'Other F', 'non_client': 'Other M'}
    colors = [F_COLOR, M_COLOR, F_COLOR_LIGHT, M_COLOR_LIGHT]
    width = 0.6

    si = sc.findfirst(sw_df.index, start_year)
    ei = sc.findfirst(sw_df.index, end_year)

    vals_all = []
    for group in groups:
        vals = np.array([
            sw_df[f'new_infections_{group}_{disease}'][si:ei].mean(),
            sw_df[f'new_transmissions_{group}_{disease}'][si:ei].mean(),
        ])
        vals_all.append(vals)
    totals = sum(vals_all)

    x = np.array([0.5, 1.5])
    bottom = np.zeros(2)
    for g, ((group, glabel), vals) in enumerate(zip(groups.items(), vals_all)):
        pct = vals / totals * 100
        p = ax.barh(x, pct, width, label=glabel, left=bottom, color=colors[g])
        if pct[0] > 5:
            ax.bar_label(p, labels=[f'{pct[0]:.0f}%', f'{pct[1]:.0f}%'],
                         label_type='center', fontsize=11)
        bottom += pct

    ax.set_title(f'{disease.upper()} infections\n{start_year}\u2013{end_year} average')
    ax.set_xlim(0, 100)
    ax.set_xlabel('%')
    ax.set_yticks(x, ['Acquired', 'Transmitted'])


def plot_hiv_prev_by_age(epi_df, ax):
    """Panel F: HIV prevalence by age and sex + ZIMPHIA"""
    thisdf = epi_df.loc[(epi_df.disease == 'hiv') &
                        (epi_df.age != '0-15') &
                        (epi_df.age != '65+')].copy()
    thisdf['prevalence'] *= 100
    sns.barplot(data=thisdf, x='age', y='prevalence', hue='sex',
                ax=ax, palette=[F_COLOR, M_COLOR], alpha=0.85)

    # Overlay ZIMPHIA data as markers
    age_groups = thisdf.age.unique()
    bar_width = 0.4
    first = True
    for i, age in enumerate(age_groups):
        if age in ZIMPHIA_HIV:
            f_val, m_val = ZIMPHIA_HIV[age]
            kw = dict(marker='D', s=60, zorder=5, edgecolors='k', linewidths=0.5)
            ax.scatter(i - bar_width/2, f_val * 100, color=F_COLOR,
                       label='ZIMPHIA 2016' if first else None, **kw)
            ax.scatter(i + bar_width/2, m_val * 100, color=M_COLOR, **kw)
            first = False

    ax.set_title('HIV prevalence by age')
    ax.set_ylabel('Prevalence (%)')
    ax.set_xlabel('Age group')
    ax.set_ylim(bottom=0)
    ax.legend(frameon=False, fontsize=11)


if __name__ == '__main__':

    epi_df = sc.loadobj('results/epi_df.df')
    sw_df = sc.loadobj('results/sw_df.df')
    cs = sc.loadobj('results/zimbabwe_calib_stats_all.df')

    set_font(size=18)
    fig = pl.figure(figsize=(22, 14))
    gs = GridSpec(2, 3, left=0.06, right=0.98, bottom=0.06, top=0.93,
                  wspace=0.28, hspace=0.38)

    ax = fig.add_subplot(gs[0, 0])
    plot_active_syph_by_hiv(cs, ax)

    ax = fig.add_subplot(gs[0, 1])
    plot_transmission_by_stage_bars(cs, ax)

    ax = fig.add_subplot(gs[0, 2])
    plot_active_prev_by_age(epi_df, ax)

    ax = fig.add_subplot(gs[1, 0])
    plot_infections_by_sw_pct(sw_df, disease='syph', ax=ax)

    ax = fig.add_subplot(gs[1, 1])
    plot_infections_by_sw_pct(sw_df, disease='hiv', ax=ax)

    ax = fig.add_subplot(gs[1, 2])
    plot_hiv_prev_by_age(epi_df, ax)

    for i, ax in enumerate(fig.axes):
        ax.text(-0.10, 1.06, chr(65 + i), transform=ax.transAxes,
                fontsize=28, fontweight='bold', va='top')

    pl.savefig('figures/fig1_syph_hiv_epi.png', dpi=200, bbox_inches='tight')
    print('Saved figures/fig1_syph_hiv_epi.png')
