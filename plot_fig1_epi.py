"""
Plot syphilis and HIV epidemiology in Zimbabwe: Fig 1 in manuscript.

Row 1 (2 panels):
    A: Active syphilis prevalence by age and sex work status + ZIMPHIA
    B: HIV prevalence by age (4 mini-panels: pre-ART, 2016, 2020, 2025)

Row 2 (3 panels):
    C: Syphilis prevalence ratio HIV+/HIV- over time
    D: Infections by sex work status (syphilis + HIV combined)
    E: Syphilis transmission by disease stage
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
# ZIMPHIA 2015-2016 HIV prevalence by age/sex (from ZIMPHIA summary sheet Fig 1)
# Format: age_group: (f_point, m_point, f_lo, f_hi, m_lo, m_hi) — 95% CIs from error bars
ZIMPHIA_HIV_2016 = {
    '15-20': (0.040, 0.025, 0.028, 0.052, 0.015, 0.035),
    '20-25': (0.077, 0.034, 0.059, 0.095, 0.020, 0.048),
    '25-30': (0.137, 0.077, 0.113, 0.161, 0.055, 0.099),
    '30-35': (0.207, 0.148, 0.175, 0.239, 0.118, 0.178),
    '35-50': (0.233, 0.207, 0.198, 0.268, 0.173, 0.241),  # avg of 35-39, 40-44, 45-49
    '50-65': (0.130, 0.142, 0.091, 0.169, 0.098, 0.186),  # avg of 50-54, 55-59, 60-64
}
# ZIMPHIA 2020 HIV prevalence by age/sex (Summary Sheet, Dec 2020)
# Format: age_group: (f_point, m_point, f_lo, f_hi, m_lo, m_hi) — 95% CIs from error bars
ZIMPHIA_HIV_2020 = {
    '15-20': (0.038, 0.021, 0.026, 0.050, 0.012, 0.030),
    '20-25': (0.064, 0.028, 0.048, 0.080, 0.016, 0.040),
    '25-30': (0.106, 0.040, 0.084, 0.128, 0.025, 0.055),
    '30-35': (0.184, 0.093, 0.154, 0.214, 0.070, 0.116),
    '35-50': (0.295, 0.208, 0.254, 0.336, 0.169, 0.247),  # avg of 35-39, 40-44, 45-49
    '50-65': (0.250, 0.251, 0.209, 0.291, 0.200, 0.302),  # avg of 50-54, 55-59, 60-64
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
    ax.set_ylabel('Prevalence (%)')
    ax.set_ylim(bottom=0)
    ax.legend(frameon=False, fontsize=14, ncol=2)


def plot_transmission_by_stage_bars(cs, ax):
    """Sexual transmission by stage (left) + MTC stacked by mother's stage and birth outcome (right)"""
    from matplotlib.patches import Patch

    # --- Sexual transmission by stage ---
    sex_stages = ['primary', 'secondary', 'early']
    sex_labels = ['Primary', 'Secondary', 'Early\nlatent']
    sex_colors = ['#e41a1c', '#ff7f00', '#984ea3']

    sex_total = 0
    sex_vals = {}
    for stage in sex_stages + ['late', 'tertiary']:
        col = f'transmission_by_stage.new_sex_{stage}'
        if col in cs.columns.get_level_values(0):
            v = cs[(col, '50%')].sum()
            sex_vals[stage] = v
            sex_total += v
    sex_pcts = [sex_vals.get(s, 0) / sex_total * 100 if sex_total > 0 else 0 for s in sex_stages]

    # --- MTC: stacked bars by mother's stage, segments = birth outcome ---
    mtc_stages = ['primary', 'secondary', 'early', 'late']
    mtc_stage_labels = ['Primary', 'Secondary', 'Early\nlatent', 'Late\nlatent']
    outcomes = ['death', 'congenital']
    outcome_colors = ['#555555', '#cb181d']
    outcome_labels = ['Adverse', 'Congenital']

    # Get MTC data
    mtc_total = 0
    mtc_data = {}
    for stage in mtc_stages:
        mtc_data[stage] = {}
        for outcome in outcomes + ['normal']:
            col = f'transmission_by_stage.new_mtc_{stage}_{outcome}'
            v = cs[(col, '50%')].sum() if col in cs.columns.get_level_values(0) else 0
            mtc_data[stage][outcome] = v
            mtc_total += v

    # --- Plot ---
    gap = 0.8
    x_sex = np.arange(len(sex_stages))
    x_mtc = np.arange(len(mtc_stages)) + len(sex_stages) + gap
    width = 0.7

    # Sexual bars
    bars_sex = ax.bar(x_sex, sex_pcts, width, color=sex_colors, alpha=0.85, edgecolor='white', linewidth=0.5)
    for bar, pct in zip(bars_sex, sex_pcts):
        if pct > 3:
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 1,
                    f'{pct:.0f}%', ha='center', va='bottom', fontsize=16, fontweight='bold', color='black')

    # MTC stacked bars — normalize within MTC so stages sum to 100%
    mtc_total_by_outcome = {}
    for outcome in outcomes:
        mtc_total_by_outcome[outcome] = sum(mtc_data[stage][outcome] for stage in mtc_stages)
    mtc_outcomes_total = sum(mtc_total_by_outcome.values())

    bottom = np.zeros(len(mtc_stages))
    for oi, (outcome, oc_color, oc_label) in enumerate(zip(outcomes, outcome_colors, outcome_labels)):
        vals = [mtc_data[stage][outcome] / mtc_outcomes_total * 100 if mtc_outcomes_total > 0 else 0 for stage in mtc_stages]
        ax.bar(x_mtc, vals, width, bottom=bottom, color=oc_color, alpha=0.85,
               edgecolor='white', linewidth=0.5, label=oc_label)
        bottom += np.array(vals)

    # Total label on top of each MTC bar
    for i, stage in enumerate(mtc_stages):
        total_pct = sum(mtc_data[stage][o] / mtc_outcomes_total * 100 for o in outcomes) if mtc_outcomes_total > 0 else 0
        if total_pct > 2:
            ax.text(x_mtc[i], bottom[i] + 0.5, f'{total_pct:.0f}%',
                    ha='center', va='bottom', fontsize=16, fontweight='bold', color='black')

    all_x = np.concatenate([x_sex, x_mtc])
    all_labels = sex_labels + mtc_stage_labels
    ax.set_xticks(all_x)
    ax.set_xticklabels(all_labels)
    ax.set_ylabel('Share of transmissions (%)')
    ax.set_title('(E) Syphilis transmission')
    ax.set_ylim(0, max(max(sex_pcts), max(bottom)) * 1.25)

    # Group labels at very top of panel (axes coords)
    ax.text(0.22, 0.98, 'Sexual', transform=ax.transAxes, ha='center', va='top', fontsize=18, fontweight='bold')
    ax.text(0.75, 0.98, 'Maternal', transform=ax.transAxes, ha='center', va='top', fontsize=18, fontweight='bold')

    # Vertical separator
    sep_x = (x_sex[-1] + x_mtc[0]) / 2
    ax.axvline(x=sep_x, color='grey', linewidth=2, linestyle='-', alpha=0.5)

    # Legend with annual counts
    recent = cs.index[(cs.index >= 2015) & (cs.index <= 2024)]
    cols = cs.columns.get_level_values(0)
    cong_col = 'syph.new_congenital'
    annual_cong = cs.loc[recent, (cong_col, '50%')].mean() if cong_col in cols else 0
    annual_death = sum(
        cs.loc[recent, (f'transmission_by_stage.new_mtc_{s}_death', '50%')].mean()
        for s in mtc_stages
        if f'transmission_by_stage.new_mtc_{s}_death' in cols
    )
    legend_labels = [
        f'Deaths ({annual_death:,.0f}/yr)',
        f'Congenital ({annual_cong:,.0f}/yr)',
    ]
    ax.legend(handles=[Patch(facecolor=c, alpha=0.85) for c in outcome_colors],
              labels=legend_labels, frameon=False, fontsize=14,
              loc='upper left', bbox_to_anchor=(0.55, 0.85))


def plot_infections_by_sw_pct_combined(cs, ax, start_year=2000, end_year=2019):
    """Infections by sex work status — syphilis and HIV side by side as vertical stacked bars"""
    groups = {'fsw': 'Women in transactional sex', 'client': 'Men in transactional sex',
               'non_fsw': 'Other girls/women', 'non_client': 'Other boys/men'}
    colors = [F_COLOR, M_COLOR, F_COLOR_LIGHT, M_COLOR_LIGHT]
    width = 0.7

    years = cs.index[(cs.index >= start_year) & (cs.index <= end_year)]

    bar_labels = ['Syph\nacquired', 'Syph\ntransmitted', 'HIV\nacquired', 'HIV\ntransmitted']
    x = np.arange(4)

    # Gather data: [syph_acquired, syph_transmitted, hiv_acquired, hiv_transmitted]
    all_vals = []
    for disease in ['syph', 'hiv']:
        for metric in ['infections', 'transmissions']:
            vals = []
            for group in groups:
                col = f'sw_stats.new_{metric}_{group}_{disease}'
                vals.append(cs.loc[years, (col, '50%')].mean())
            all_vals.append(np.array(vals))

    # Normalize each bar to 100%
    totals = [v.sum() for v in all_vals]

    bottom = np.zeros(4)
    for g, (group, glabel) in enumerate(groups.items()):
        pcts = np.array([all_vals[i][g] / totals[i] * 100 if totals[i] > 0 else 0 for i in range(4)])
        ax.bar(x, pcts, width, label=glabel, bottom=bottom, color=colors[g])
        for j in range(4):
            if pcts[j] > 5:
                ax.text(x[j], bottom[j] + pcts[j]/2, f'{pcts[j]:.0f}%',
                        ha='center', va='center', fontsize=16, color='black', fontweight='bold')
        bottom += pcts

    ax.set_title(f'(D) Infections by sex work\n{start_year}\u2013{end_year} avg')
    ax.set_ylim(0, 120)
    ax.set_xticks(x)
    ax.set_xticklabels(bar_labels)

    # Vertical separator between syphilis and HIV bars
    ax.axvline(x=1.5, color='grey', linewidth=1.2, linestyle='-', alpha=0.4)

    ax.legend(loc='upper left', fontsize=14, frameon=False,
              ncol=2, columnspacing=1.0, bbox_to_anchor=(0.0, 1.0))


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
    ax.legend(frameon=False, fontsize=14)
    ax.set_ylabel('Prevalence ratio')
    ax.set_title('(C) Syphilis prevalence\nHIV+ / HIV\u2212')
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

    ymax = 50  # Shared y-axis

    for idx, (year, label, zimphia, marker) in enumerate(year_specs):
        ax = fig.add_subplot(inner_gs[idx])

        if idx == 0:
            ax.set_title('(B) HIV prevalence by age and year')

        if year not in cs.index:
            ax.text(0.5, 0.5, f'{year}: no data', transform=ax.transAxes, ha='center')
            continue

        f_vals, m_vals = [], []
        f_errs, m_errs = [[], []], [[], []]
        for lo, hi in age_bins:
            f_col = f'epi_ts.hiv_prev_f_{lo}_{hi}'
            m_col = f'epi_ts.hiv_prev_m_{lo}_{hi}'
            cols = cs.columns.get_level_values(0)
            f_med = cs.loc[year, (f_col, '50%')] * 100 if f_col in cols else 0
            m_med = cs.loc[year, (m_col, '50%')] * 100 if m_col in cols else 0
            f_lo = cs.loc[year, (f_col, '10%')] * 100 if f_col in cols else 0
            f_hi = cs.loc[year, (f_col, '90%')] * 100 if f_col in cols else 0
            m_lo = cs.loc[year, (m_col, '10%')] * 100 if m_col in cols else 0
            m_hi = cs.loc[year, (m_col, '90%')] * 100 if m_col in cols else 0
            f_vals.append(f_med); m_vals.append(m_med)
            f_errs[0].append(f_med - f_lo); f_errs[1].append(f_hi - f_med)
            m_errs[0].append(m_med - m_lo); m_errs[1].append(m_hi - m_med)

        f_label = 'Female' if idx == 0 else None
        m_label = 'Male' if idx == 0 else None
        ax.bar(x - width/2, f_vals, width, color=F_COLOR, alpha=0.85, edgecolor='white', linewidth=0.5, label=f_label)
        ax.bar(x + width/2, m_vals, width, color=M_COLOR, alpha=0.85, edgecolor='white', linewidth=0.5, label=m_label)
        ax.errorbar(x - width/2, f_vals, yerr=f_errs, fmt='none', ecolor='grey', capsize=2, linewidth=0.7, alpha=0.5)
        ax.errorbar(x + width/2, m_vals, yerr=m_errs, fmt='none', ecolor='grey', capsize=2, linewidth=0.7, alpha=0.5)

        # ZIMPHIA overlay with CIs
        if zimphia is not None:
            kw = dict(s=50, zorder=5, edgecolors='k', linewidths=0.5, marker=marker)
            first_zimphia = True
            for j, al in enumerate(age_labels):
                if al in zimphia:
                    vals = zimphia[al]
                    f_val, m_val = vals[0], vals[1]
                    zlabel = f'ZIMPHIA {label}' if first_zimphia else None
                    ax.scatter(x[j] - width/2, f_val * 100, color=F_COLOR, label=zlabel, **kw)
                    ax.scatter(x[j] + width/2, m_val * 100, color=M_COLOR, **kw)
                    first_zimphia = False
                    if len(vals) == 6:
                        f_lo, f_hi, m_lo, m_hi = vals[2], vals[3], vals[4], vals[5]
                        ax.errorbar(x[j] - width/2, f_val * 100,
                                    yerr=[[(f_val - f_lo) * 100], [(f_hi - f_val) * 100]],
                                    fmt='none', ecolor=F_COLOR, capsize=2, linewidth=1, alpha=0.7,
                                    elinewidth=1.5, capthick=1, zorder=4)
                        ax.errorbar(x[j] + width/2, m_val * 100,
                                    yerr=[[(m_val - m_lo) * 100], [(m_hi - m_val) * 100]],
                                    fmt='none', ecolor=M_COLOR, capsize=2, linewidth=1, alpha=0.7,
                                    elinewidth=1.5, capthick=1, zorder=4)

        ax.set_ylim(0, ymax)
        ax.set_xlim(-0.6, len(age_bins) - 0.4)

        if idx == 0:
            ax.legend(frameon=False, fontsize=14, loc='upper left', ncol=2)
        elif idx in [1, 2]:
            ax.legend(frameon=False, fontsize=14, loc='upper left')

        # Year label inside panel
        ax.text(0.98, 0.85, label, transform=ax.transAxes, ha='right', va='top',
                fontsize=16, fontweight='bold')

        # Only show x labels on bottom panel
        if idx == len(year_specs) - 1:
            ax.set_xticks(x)
            ax.set_xticklabels(age_labels)
        else:
            ax.set_xticks(x)
            ax.set_xticklabels([])

        # Only show y label on middle panels
        if idx == 1:
            ax.set_ylabel('HIV prevalence (%)')


if __name__ == '__main__':

    sw_prev_df = sc.loadobj('results/sw_prev_df.df')
    cs = sc.loadobj('results/zimbabwe_calib_stats_all.df')

    set_font(size=20)
    from matplotlib.gridspec import GridSpecFromSubplotSpec

    fig = pl.figure(figsize=(22, 16))
    gs = GridSpec(2, 1, left=0.04, right=0.98, bottom=0.05, top=0.95,
                  hspace=0.22, height_ratios=[1, 1])

    # --- Row 1: Prevalence by age (2 panels) ---
    gs_row1 = GridSpecFromSubplotSpec(1, 2, subplot_spec=gs[0], wspace=0.15,
                                      width_ratios=[1.9, 1.3])

    ax_a = fig.add_subplot(gs_row1[0, 0])
    plot_prev_by_sw(sw_prev_df, ax=ax_a, prev_col='active_prevalence',
                    title='(A) Active syphilis prevalence\nby age and sex work status',
                    show_zimphia=True, zimphia_by_age=ZIMPHIA_SYPH_ACTIVE_BY_AGE,
                    exclude_ages=['50-65'])

    # Panel B: HIV prev by age — 4 stacked mini-panels
    gs_b = GridSpecFromSubplotSpec(4, 1, subplot_spec=gs_row1[0, 1], hspace=0.18)
    plot_hiv_prev_by_age_lines(cs, fig, gs_b)

    # --- Row 2: 3 panels (C, D, E) ---
    gs_row2 = GridSpecFromSubplotSpec(1, 3, subplot_spec=gs[1], wspace=0.15,
                                      width_ratios=[0.7, 1, 1.5])

    ax_c = fig.add_subplot(gs_row2[0, 0])
    plot_syph_hiv_ratio_ts(cs, ax_c, start_year=2000)

    ax_d = fig.add_subplot(gs_row2[0, 1])
    plot_infections_by_sw_pct_combined(cs, ax=ax_d)

    ax_e = fig.add_subplot(gs_row2[0, 2])
    plot_transmission_by_stage_bars(cs, ax_e)

    pl.savefig('figures/fig1_syph_hiv_epi.png', dpi=200, bbox_inches='tight')
    print('Saved figures/fig1_syph_hiv_epi.png')
