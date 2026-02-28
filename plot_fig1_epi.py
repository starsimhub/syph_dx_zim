"""
Plot syphilis and HIV epidemiology in Zimbabwe: Fig 1 in manuscript.

Row 1 (3 panels):
    A: Active syphilis prevalence by age and sex work status + ZIMPHIA
    B: Syphilis transmission by disease stage
    C: Syphilis infections by sex work status (%)

Row 2 (3 panels):
    D: Syphilis prevalence ratio HIV+/HIV- over time
    E: HIV infections by sex work status (%)
    F: HIV prevalence by age + ZIMPHIA 2016 & 2020
"""

import sciris as sc
import numpy as np
import pandas as pd
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
                    show_zimphia=False, zimphia_by_age=None, exclude_ages=None):
    """Prevalence by age with 3 bars (FSW/clients, all F, all M) + aggregated 15-65 bar"""
    skip = {'0-15', '65+'} | set(exclude_ages or [])
    df = sw_prev_df[(sw_prev_df.disease == 'syph') &
                    (~sw_prev_df.age.isin(skip))]

    age_groups = [a for a in df.age.unique() if a not in skip]

    # 3 bars: SW women, all women, all men
    bar_specs = [
        ('Female', 'SW',      F_COLOR,   'Women in transactional sex'),
        ('Female', 'Overall', F_COLOR_LIGHT, 'All girls/women'),
        ('Male',   'Overall', M_COLOR_LIGHT, 'All boys/men'),
    ]
    n_bars = len(bar_specs)
    width = 0.25
    # Aggregate label reflects the actual age range included
    first_lo = age_groups[0].split('-')[0] if age_groups else '15'
    last_hi = age_groups[-1].split('-')[1] if age_groups else '65'
    all_labels = age_groups + [f'{first_lo}-{last_hi}']
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


def plot_transmission_by_stage_bars(cs, ax):
    """Grouped bars: % of sexual and MTC transmission by stage"""
    stages = ['primary', 'secondary', 'early']
    stage_labels = ['Primary', 'Secondary', 'Early\nlatent']
    colors = ['#e41a1c', '#ff7f00', '#984ea3']

    # Compute percentages for each transmission mode
    pcts_by_mode = {}
    for mode in ['sex', 'mtc']:
        total = 0
        vals = {}
        for stage in stages + ['late', 'tertiary']:
            col = f'transmission_by_stage.new_{mode}_{stage}'
            if col in cs.columns.get_level_values(0):
                v = cs[(col, '50%')].sum()
                vals[stage] = v
                total += v
        if total > 0:
            pcts_by_mode[mode] = [vals.get(s, 0) / total * 100 for s in stages]
        else:
            pcts_by_mode[mode] = [0] * len(stages)

    x = np.arange(len(stages))
    width = 0.35

    # Sexual transmission bars
    bars_sex = ax.bar(x - width/2, pcts_by_mode['sex'], width, color=colors,
                      alpha=0.85, edgecolor='white', linewidth=0.5)
    # MTC bars (lighter)
    bars_mtc = ax.bar(x + width/2, pcts_by_mode['mtc'], width, color=colors,
                      alpha=0.4, edgecolor='white', linewidth=0.5, hatch='///')

    for bar, pct in zip(bars_sex, pcts_by_mode['sex']):
        if pct > 3:
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 1,
                    f'{pct:.0f}%', ha='center', va='bottom', fontsize=11, fontweight='bold')
    for bar, pct in zip(bars_mtc, pcts_by_mode['mtc']):
        if pct > 3:
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 1,
                    f'{pct:.0f}%', ha='center', va='bottom', fontsize=11, fontweight='bold')

    ax.set_xticks(x)
    ax.set_xticklabels(stage_labels)
    ax.set_ylabel('Share of\ntransmissions (%)')
    ax.set_title('Syphilis transmission\nby disease stage')
    ax.set_ylim(0, 80)

    # Legend for transmission mode
    from matplotlib.patches import Patch
    ax.legend([Patch(facecolor='grey', alpha=0.85), Patch(facecolor='grey', alpha=0.4, hatch='///')],
              ['Sexual', 'Maternal'], frameon=False, fontsize=12, loc='upper right')


def plot_infections_by_sw_pct(sw_df, disease, ax, start_year=2000, end_year=2019):
    """Infections by sex work status as percentages"""
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

    disease_label = 'Syphilis' if disease == 'syph' else 'HIV'
    ax.set_title(f'{disease_label} infections\n{start_year}\u2013{end_year} avg')
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 2.3)
    ax.set_xlabel('%')
    ax.set_yticks(x, ['Acquired', 'Transmitted'])

    ax.legend(loc='upper right', fontsize=10, frameon=False,
              ncol=2, columnspacing=0.8, handlelength=1.0, labelspacing=0.3)


def plot_syph_hiv_ratio_ts(cs, ax, start_year=1995):
    """Ratio of syphilis prevalence in HIV+ vs HIV- over time"""
    years = cs.index[cs.index >= start_year]

    hiv_pos = cs.loc[years, ('coinfection_stats.syph_prev_has_hiv', '50%')]
    hiv_neg = cs.loc[years, ('coinfection_stats.syph_prev_no_hiv', '50%')]
    ratio_med = hiv_pos / hiv_neg

    hiv_pos_lo = cs.loc[years, ('coinfection_stats.syph_prev_has_hiv', '10%')]
    hiv_neg_hi = cs.loc[years, ('coinfection_stats.syph_prev_no_hiv', '90%')]
    hiv_pos_hi = cs.loc[years, ('coinfection_stats.syph_prev_has_hiv', '90%')]
    hiv_neg_lo = cs.loc[years, ('coinfection_stats.syph_prev_no_hiv', '10%')]
    ratio_lo = hiv_pos_lo / hiv_neg_hi
    ratio_hi = hiv_pos_hi / hiv_neg_lo

    ax.fill_between(years, ratio_lo, ratio_hi, alpha=0.15, color='#555555')
    ax.plot(years, ratio_med, color='#555555', linewidth=2, label='Model')

    # ZIMPHIA 2016 ratio: 2.9% / 0.4% = 7.25
    zimphia_ratio = ZIMPHIA_SYPH['active_hiv_pos'] / ZIMPHIA_SYPH['active_hiv_neg']
    ax.scatter([2016], [zimphia_ratio], color=HIV_POS_COLOR,
               marker='D', s=80, zorder=5, edgecolors='k', linewidths=0.5,
               label=f'ZIMPHIA 2016 ({zimphia_ratio:.1f}x)')

    ax.axhline(y=1, color='grey', linestyle='--', linewidth=0.8, alpha=0.5)
    ax.legend(frameon=False, fontsize=11)
    ax.set_ylabel('Prevalence ratio')
    ax.set_title('Syphilis prevalence\nHIV+ / HIV\u2212')
    ax.set_ylim(bottom=0)
    ax.set_xlim(start_year, 2025)


def plot_hiv_prev_by_age_lines(cs, fig, inner_gs):
    """HIV prevalence by age as 4 stacked mini-panels, one per year snapshot"""
    age_bins = [(15, 20), (20, 25), (25, 30), (30, 35), (35, 50), (50, 65)]
    age_labels = ['15-20', '20-25', '25-30', '30-35', '35-50', '50-65']
    x = np.arange(len(age_bins))
    width = 0.35

    year_specs = [
        (2005, 'Pre-ART', None, None),
        (2016, '2016', ZIMPHIA_HIV_2016, 'D'),
        (2020, '2020', ZIMPHIA_HIV_2020, 's'),
        (2025, '2025 (proj.)', None, None),
    ]

    ymax = 40  # Shared y-axis

    for idx, (year, label, zimphia, marker) in enumerate(year_specs):
        ax = fig.add_subplot(inner_gs[idx])

        if year not in cs.index:
            ax.text(0.5, 0.5, f'{year}: no data', transform=ax.transAxes, ha='center')
            continue

        f_vals = []
        m_vals = []
        for lo, hi in age_bins:
            f_col = f'epi_ts.hiv_prev_f_{lo}_{hi}'
            m_col = f'epi_ts.hiv_prev_m_{lo}_{hi}'
            f_vals.append(cs.loc[year, (f_col, '50%')] * 100 if f_col in cs.columns.get_level_values(0) else 0)
            m_vals.append(cs.loc[year, (m_col, '50%')] * 100 if m_col in cs.columns.get_level_values(0) else 0)

        ax.bar(x - width/2, f_vals, width, color=F_COLOR, alpha=0.85, edgecolor='white', linewidth=0.5)
        ax.bar(x + width/2, m_vals, width, color=M_COLOR, alpha=0.85, edgecolor='white', linewidth=0.5)

        # ZIMPHIA overlay
        if zimphia is not None:
            kw = dict(s=50, zorder=5, edgecolors='k', linewidths=0.5, marker=marker)
            for j, al in enumerate(age_labels):
                if al in zimphia:
                    f_val, m_val = zimphia[al]
                    ax.scatter(x[j] - width/2, f_val * 100, color=F_COLOR, **kw)
                    ax.scatter(x[j] + width/2, m_val * 100, color=M_COLOR, **kw)

        ax.set_ylim(0, ymax)
        ax.set_xlim(-0.6, len(age_bins) - 0.4)

        # Year label inside panel
        ax.text(0.98, 0.85, label, transform=ax.transAxes, ha='right', va='top',
                fontsize=13, fontweight='bold')

        # Only show x labels on bottom panel
        if idx == len(year_specs) - 1:
            ax.set_xticks(x)
            ax.set_xticklabels(age_labels, fontsize=11)
            ax.set_xlabel('Age group')
        else:
            ax.set_xticks(x)
            ax.set_xticklabels([])

        # Only show y label on middle panels
        if idx == 1:
            ax.set_ylabel('HIV prevalence (%)')

        ax.tick_params(labelsize=11)


if __name__ == '__main__':

    epi_df = sc.loadobj('results/epi_df.df')
    sw_df = sc.loadobj('results/sw_df.df')
    sw_prev_df = sc.loadobj('results/sw_prev_df.df')
    cs = sc.loadobj('results/zimbabwe_calib_stats_all.df')

    set_font(size=18)
    fig = pl.figure(figsize=(22, 14))
    gs = GridSpec(2, 3, left=0.06, right=0.98, bottom=0.06, top=0.93,
                  wspace=0.28, hspace=0.38)

    # --- Row 1: Syphilis ---
    ax = fig.add_subplot(gs[0, 0])
    plot_prev_by_sw(sw_prev_df, ax=ax, prev_col='active_prevalence',
                    title='Active syphilis prevalence\nby sex work status',
                    show_zimphia=True, zimphia_by_age=ZIMPHIA_SYPH_ACTIVE_BY_AGE,
                    exclude_ages=['50-65'])

    ax = fig.add_subplot(gs[0, 1])
    plot_infections_by_sw_pct(sw_df, disease='syph', ax=ax)

    ax = fig.add_subplot(gs[0, 2])
    plot_transmission_by_stage_bars(cs, ax)

    # --- Row 2: HIV + coinfection ---
    # Panel D: 4 stacked mini-panels for HIV by age at different years
    from matplotlib.gridspec import GridSpecFromSubplotSpec
    inner_gs = GridSpecFromSubplotSpec(4, 1, subplot_spec=gs[1, 0], hspace=0.08)
    plot_hiv_prev_by_age_lines(cs, fig, inner_gs)

    ax = fig.add_subplot(gs[1, 1])
    plot_infections_by_sw_pct(sw_df, disease='hiv', ax=ax)

    ax = fig.add_subplot(gs[1, 2])
    plot_syph_hiv_ratio_ts(cs, ax)

    # Label panels: A-C for row 1, D for the stacked mini-panels, E-F for rest of row 2
    main_axes = [fig.axes[0], fig.axes[1], fig.axes[2]]  # A, B, C
    # The inner mini-panels are axes 3-6; label only the first one as D
    main_axes.append(fig.axes[3])  # D (top mini-panel)
    # E and F are the last two main axes
    main_axes.append(fig.axes[7])  # E (HIV infections)
    main_axes.append(fig.axes[8])  # F (ratio)
    for i, ax in enumerate(main_axes):
        ax.text(-0.10, 1.06, chr(65 + i), transform=ax.transAxes,
                fontsize=28, fontweight='bold', va='top')

    pl.savefig('figures/fig1_syph_hiv_epi.png', dpi=200, bbox_inches='tight')
    print('Saved figures/fig1_syph_hiv_epi.png')
