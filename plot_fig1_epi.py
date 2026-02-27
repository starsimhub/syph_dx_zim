"""
Plot syphilis and HIV epidemiology in Zimbabwe: Fig 1 in manuscript.

Row 1 (3 panels, equal width):
    A: Active syphilis prevalence by age and sex work status + ZIMPHIA
    B: All syphilis prevalence by age and sex work status + WHO
    C: Ever exposed to syphilis (time series by sex/FSW)

Row 2 (4 panels: 3 narrow + 1 wider):
    D: Active syphilis prevalence by HIV status (bars + ZIMPHIA)
    E: Syphilis transmission by disease stage
    F: Syphilis infections by sex work status (%)
    G: HIV prevalence by age + ZIMPHIA 2016 & 2020
"""

import sciris as sc
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as pl
from matplotlib.gridspec import GridSpec
from utils import set_font


# ZIMPHIA 2015-2016 data — overall
ZIMPHIA_SYPH = dict(
    active_hiv_pos=0.029, active_hiv_neg=0.004,
    active_f=0.010, active_m=0.006,
    ever_exposed_f=0.023, ever_exposed_m=0.019,  # Overall 15-49 (HIV+ and HIV- combined)
)
# ZIMPHIA 2015-2016 syphilis by age/sex (from published ZIMPHIA report)
# "active" = treponemal+ AND RPR+ (includes latent → maps to model "all infected")
# "ever" = treponemal+ only (maps to model "ever_exposed")
ZIMPHIA_SYPH_ACTIVE_BY_AGE = {
    # age_group: (female, male)
    '15-20': (0.004, 0.002),
    '20-25': (0.013, 0.004),
    '25-30': (0.011, 0.004),
    '30-35': (0.010, 0.008),
    '35-50': (0.008, 0.009),   # avg of 35-39, 40-44, 45-49
    '50-65': (0.014, 0.008),   # avg of 50-54, 55-59, 60-64
}
ZIMPHIA_SYPH_EVER_BY_AGE = {
    # age_group: (female, male)
    '15-20': (0.008, 0.003),
    '20-25': (0.021, 0.010),
    '25-30': (0.024, 0.013),
    '30-35': (0.021, 0.018),
    '35-50': (0.034, 0.036),   # avg of 35-39, 40-44, 45-49
    '50-65': (0.086, 0.077),   # avg of 50-54, 55-59, 60-64
}
# WHO GHO: syphilis seroprevalence among FSW in Zimbabwe (~30%, treponemal = ever-exposed)
WHO_FSW_SYPH_PREV = 0.30
# ZIMPHIA 2015-2016 HIV prevalence by age/sex (from published ZIMPHIA report)
ZIMPHIA_HIV_2016 = {
    # age_group: (female, male)
    '15-20': (0.040, 0.025),
    '20-25': (0.077, 0.034),
    '25-30': (0.137, 0.077),
    '30-35': (0.207, 0.148),
    '35-50': (0.233, 0.207),
    '50-65': (0.130, 0.142),
}
# ZIMPHIA 2020 HIV prevalence by age/sex (Summary Sheet, Dec 2020)
# Mapped to model age bins by averaging 5-year ZIMPHIA groups
ZIMPHIA_HIV_2020 = {
    # age_group: (female, male)
    '15-20': (0.038, 0.021),
    '20-25': (0.064, 0.028),
    '25-30': (0.106, 0.040),
    '30-35': (0.184, 0.093),
    '35-50': (0.295, 0.208),   # avg of 35-39, 40-44, 45-49
    '50-65': (0.250, 0.251),   # avg of 50-54, 55-59, 60-64
}

# Colors — sex-disaggregated
F_COLOR = '#d46e9c'
M_COLOR = '#4a90d9'
F_COLOR_LIGHT = '#f0a3c4'
M_COLOR_LIGHT = '#a3c4e8'

# Colors — HIV-disaggregated
HIV_POS_COLOR = '#e41a1c'
HIV_NEG_COLOR = '#377eb8'


def plot_prev_by_sw(sw_prev_df, ax, prev_col='active_prevalence', title=None,
                    show_zimphia=False, zimphia_by_age=None):
    """Prevalence by age with 3 bars (FSW/clients, all F, all M) + aggregated 15-65 bar"""
    df = sw_prev_df[(sw_prev_df.disease == 'syph') &
                    (sw_prev_df.age != '0-15') &
                    (sw_prev_df.age != '65+')]

    age_groups = [a for a in df.age.unique() if a not in ('0-15', '65+')]

    # 3 bars: SW women, all women, all men
    bar_specs = [
        ('Female', 'SW',      F_COLOR,   'Women in transactional sex'),
        ('Female', 'Overall', F_COLOR_LIGHT, 'All girls/women'),
        ('Male',   'Overall', M_COLOR_LIGHT, 'All boys/men'),
    ]
    n_bars = len(bar_specs)
    width = 0.25
    all_labels = age_groups + ['15-65']
    x = np.arange(len(all_labels))
    agg_x = x[-1]

    for i, (sex, sw_group, color, label) in enumerate(bar_specs):
        vals = []
        errs = []
        for age in age_groups:
            sub = df[(df.age == age) & (df.sex == sex) & (df.sw_group == sw_group)]
            vals.append(sub[prev_col].mean() * 100 if len(sub) > 0 else 0)
            errs.append(sub[prev_col].std() * 100 if len(sub) > 1 else 0)
        # Aggregate 15-65
        agg_sub = df[(df.sex == sex) & (df.sw_group == sw_group)]
        agg_val = agg_sub.groupby('par_idx')[prev_col].mean() * 100 if len(agg_sub) > 0 else pd.Series([0])
        vals.append(agg_val.mean())
        errs.append(agg_val.std() if len(agg_val) > 1 else 0)

        offset = (i - (n_bars - 1) / 2) * width
        ax.bar(x + offset, vals, width, color=color, alpha=0.85,
               edgecolor='white', linewidth=0.5, label=label)
        ax.errorbar(x + offset, vals, yerr=errs, fmt='none', ecolor='grey',
                    capsize=2, linewidth=0.8, alpha=0.7)

    # Vertical separator before aggregate bar
    ax.axvline(x=agg_x - 0.6, color='grey', linewidth=0.5, linestyle='-', alpha=0.3)

    # ZIMPHIA data overlays
    kw = dict(marker='D', s=50, zorder=5, edgecolors='k', linewidths=0.5)
    f_offset = (1 - (n_bars - 1) / 2) * width  # All F bar position
    m_offset = (2 - (n_bars - 1) / 2) * width  # All M bar position

    if zimphia_by_age is not None:
        first = True
        for j, age in enumerate(age_groups):
            if age in zimphia_by_age:
                f_val, m_val = zimphia_by_age[age]
                ax.scatter(x[j] + f_offset, f_val * 100, color=F_COLOR_LIGHT,
                           label='ZIMPHIA 2016' if first else None, **kw)
                ax.scatter(x[j] + m_offset, m_val * 100, color=M_COLOR_LIGHT, **kw)
                first = False

    if show_zimphia:
        ax.scatter(agg_x + f_offset, ZIMPHIA_SYPH['active_f'] * 100, color=F_COLOR_LIGHT,
                   label='ZIMPHIA 2016' if zimphia_by_age is None else None, **kw)
        ax.scatter(agg_x + m_offset, ZIMPHIA_SYPH['active_m'] * 100, color=M_COLOR_LIGHT, **kw)

    if title is None:
        title = 'Syphilis prevalence\nby sex work status'
    ax.set_title(title)
    ax.set_xticks(x)
    ax.set_xticklabels(all_labels)
    ax.set_xlabel('Age group')
    ax.set_ylabel('Prevalence (%)')
    ax.set_ylim(bottom=0)
    ax.legend(frameon=False, fontsize=9, ncol=2)


def plot_ever_exposed_ts(cs, ax, start_year=1995):
    """Panel C: Ever exposed to syphilis over time, by sex and FSW status"""
    years = cs.index[cs.index >= start_year]

    for col, color, label in [
        ('epi_ts.syph_ever_exposed_fsw', F_COLOR, 'FSW'),
        ('epi_ts.syph_ever_exposed_f', F_COLOR_LIGHT, 'All girls/women'),
        ('epi_ts.syph_ever_exposed_m', M_COLOR_LIGHT, 'All boys/men'),
    ]:
        med = cs.loc[years, (col, '50%')] * 100
        lo = cs.loc[years, (col, '10%')] * 100
        hi = cs.loc[years, (col, '90%')] * 100
        ax.fill_between(years, lo, hi, alpha=0.15, color=color)
        ax.plot(years, med, color=color, label=label, linewidth=2)

    # ZIMPHIA 2016 ever-exposed data points
    kw = dict(marker='D', s=80, zorder=5, edgecolors='k', linewidths=0.5)
    ax.scatter([2016], [ZIMPHIA_SYPH['ever_exposed_f'] * 100], color=F_COLOR_LIGHT,
               label='ZIMPHIA 2016', **kw)
    ax.scatter([2016], [ZIMPHIA_SYPH['ever_exposed_m'] * 100], color=M_COLOR_LIGHT, **kw)

    # WHO GHO FSW syphilis seroprevalence (~30%, treponemal = ever-exposed)
    ax.scatter([2016], [WHO_FSW_SYPH_PREV * 100], color=F_COLOR,
               marker='D', s=80, zorder=5, edgecolors='k', linewidths=0.5, label='WHO GHO FSW')

    ax.legend(frameon=False, fontsize=11, loc='upper left')
    ax.set_ylabel('Prevalence (%)')
    ax.set_title('Ever exposed to syphilis')
    ax.set_ylim(bottom=0)
    ax.set_xlim(start_year, 2025)


def plot_syph_prev_by_hiv_bars(cs, ax, year=2016):
    """Panel D: Active syphilis prevalence by HIV status — bars with error bars + ZIMPHIA"""
    x = np.array([0, 1])
    labels = ['HIV+', 'HIV\u2212']
    colors = [HIV_POS_COLOR, HIV_NEG_COLOR]
    cols = ['active_coinfection_stats.syph_prev_has_hiv', 'active_coinfection_stats.syph_prev_no_hiv']

    meds = []
    los = []
    his = []
    for col in cols:
        med = cs.loc[year, (col, '50%')] * 100
        lo = cs.loc[year, (col, '10%')] * 100
        hi = cs.loc[year, (col, '90%')] * 100
        meds.append(med)
        los.append(med - lo)
        his.append(hi - med)

    bars = ax.bar(x, meds, color=colors, alpha=0.85, edgecolor='white', linewidth=0.5, width=0.6)
    ax.errorbar(x, meds, yerr=[los, his], fmt='none', ecolor='k', capsize=5, linewidth=1.5)

    # ZIMPHIA data as diamonds
    kw = dict(marker='D', s=80, zorder=5, edgecolors='k', linewidths=0.5)
    ax.scatter(0, ZIMPHIA_SYPH['active_hiv_pos'] * 100, color=HIV_POS_COLOR,
               label='ZIMPHIA 2016', **kw)
    ax.scatter(1, ZIMPHIA_SYPH['active_hiv_neg'] * 100, color=HIV_NEG_COLOR, **kw)

    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel('Active syphilis\nprevalence (%)')
    ax.set_title('Syphilis by\nHIV status')
    ax.set_ylim(bottom=0)
    ax.legend(frameon=False, fontsize=11, loc='upper right')


def plot_transmission_by_stage_bars(cs, ax):
    """Panel E: 3 bars showing % of sexual transmission by stage"""
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
                f'{pct:.0f}%', ha='center', va='bottom', fontsize=13, fontweight='bold')

    ax.set_ylabel('Share of sexual\ntransmissions (%)')
    ax.set_title('Transmission\nby stage')
    ax.set_ylim(0, 80)
    ax.set_xlabel('')


def plot_infections_by_sw_pct(sw_df, disease, ax, start_year=2000, end_year=2019):
    """Panel F: Infections by sex work status as percentages"""
    groups = {'fsw': 'Women in transactional sex', 'client': 'Men in transactional sex',
               'non_fsw': 'Other girls/women', 'non_client': 'Other boys/men'}
    colors = [F_COLOR, M_COLOR, F_COLOR_LIGHT, M_COLOR_LIGHT]
    width = 0.85

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
                         label_type='center', fontsize=10)
        bottom += pct

    ax.set_title(f'Infections\n{start_year}\u2013{end_year} avg')
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 2.3)
    ax.set_xlabel('%')
    ax.set_yticks(x, ['Acquired', 'Transmitted'])

    ax.legend(loc='upper right', fontsize=10, frameon=False,
              ncol=2, columnspacing=0.8, handlelength=1.0, labelspacing=0.3)


def plot_hiv_prev_by_age(epi_df, ax):
    """Panel G: HIV prevalence by age and sex + ZIMPHIA 2016 & 2020"""
    thisdf = epi_df.loc[(epi_df.disease == 'hiv') &
                        (epi_df.age != '0-15') &
                        (epi_df.age != '65+')].copy()
    thisdf['prevalence'] *= 100
    sns.barplot(data=thisdf, x='age', y='prevalence', hue='sex',
                ax=ax, palette=[F_COLOR, M_COLOR], alpha=0.85)

    # Overlay ZIMPHIA data as markers
    age_groups = thisdf.age.unique()
    bar_width = 0.4
    kw_2016 = dict(marker='D', s=60, zorder=5, edgecolors='k', linewidths=0.5)
    kw_2020 = dict(marker='s', s=60, zorder=5, edgecolors='k', linewidths=0.5)
    first_2016 = True
    first_2020 = True
    for i, age in enumerate(age_groups):
        if age in ZIMPHIA_HIV_2016:
            f_val, m_val = ZIMPHIA_HIV_2016[age]
            ax.scatter(i - bar_width/2, f_val * 100, color=F_COLOR,
                       label='ZIMPHIA 2016' if first_2016 else None, **kw_2016)
            ax.scatter(i + bar_width/2, m_val * 100, color=M_COLOR, **kw_2016)
            first_2016 = False
        if age in ZIMPHIA_HIV_2020:
            f_val, m_val = ZIMPHIA_HIV_2020[age]
            ax.scatter(i - bar_width/2, f_val * 100, color=F_COLOR,
                       label='ZIMPHIA 2020' if first_2020 else None, **kw_2020)
            ax.scatter(i + bar_width/2, m_val * 100, color=M_COLOR, **kw_2020)
            first_2020 = False

    ax.set_title('HIV prevalence by age')
    ax.set_ylabel('Prevalence (%)')
    ax.set_xlabel('Age group')
    ax.set_ylim(bottom=0)
    ax.legend(frameon=False, fontsize=10)


if __name__ == '__main__':

    epi_df = sc.loadobj('results/epi_df.df')
    sw_df = sc.loadobj('results/sw_df.df')
    sw_prev_df = sc.loadobj('results/sw_prev_df.df')
    cs = sc.loadobj('results/zimbabwe_calib_stats_all.df')

    set_font(size=18)
    fig = pl.figure(figsize=(24, 14))

    # Row 1: 3 equal-width syphilis by-age panels
    gs1 = GridSpec(1, 3, left=0.05, right=0.98, bottom=0.55, top=0.93, wspace=0.28)

    # Row 2: 3 narrow panels + 1 wider panel
    gs2 = GridSpec(1, 4, left=0.05, right=0.98, bottom=0.06, top=0.44,
                   wspace=0.32, width_ratios=[1, 1, 1, 1.5])

    # --- Row 1 ---
    # Panel A: active syphilis (primary+secondary) + ZIMPHIA "active" by age
    ax = fig.add_subplot(gs1[0, 0])
    plot_prev_by_sw(sw_prev_df, ax=ax, prev_col='active_prevalence',
                    title='Active syphilis prevalence\nby sex work status',
                    show_zimphia=True, zimphia_by_age=ZIMPHIA_SYPH_ACTIVE_BY_AGE)

    # Panel B: all syphilis (any stage including latent)
    ax = fig.add_subplot(gs1[0, 1])
    plot_prev_by_sw(sw_prev_df, ax=ax, prev_col='prevalence',
                    title='All syphilis prevalence\nby sex work status')

    ax = fig.add_subplot(gs1[0, 2])
    plot_ever_exposed_ts(cs, ax)

    # --- Row 2 ---
    ax = fig.add_subplot(gs2[0, 0])
    plot_syph_prev_by_hiv_bars(cs, ax)

    ax = fig.add_subplot(gs2[0, 1])
    plot_transmission_by_stage_bars(cs, ax)

    ax = fig.add_subplot(gs2[0, 2])
    plot_infections_by_sw_pct(sw_df, disease='syph', ax=ax)

    ax = fig.add_subplot(gs2[0, 3])
    plot_hiv_prev_by_age(epi_df, ax)

    for i, ax in enumerate(fig.axes):
        ax.text(-0.10, 1.06, chr(65 + i), transform=ax.transAxes,
                fontsize=28, fontweight='bold', va='top')

    pl.savefig('figures/fig1_syph_hiv_epi.png', dpi=200, bbox_inches='tight')
    print('Saved figures/fig1_syph_hiv_epi.png')
