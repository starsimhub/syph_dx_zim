"""
Step-down prevalence figure: seroprevalence → active → primary → visible ulcerative.

Shows the dramatic reduction from seroprevalence (what an RDT detects) down to
visible ulcerative prevalence (what presents as GUD), highlighting that improving
GUD diagnostic sensitivity addresses only a small part of the bottleneck.

Structure matches Fig 1A: age groups on x-axis, grouped bars for SW women,
all women, all men.
"""

import numpy as np
import sciris as sc
import pylab as pl

RESULTS_DIR = 'results'
FIGURES_DIR = 'figures'

# Duration-based estimate: primary fraction of active (primary + secondary)
# Primary ~6 weeks, secondary ~14.4 weeks (3.6 months)
FRAC_PRIMARY = 6 / (6 + 14.4)  # ~0.294

# Visibility of primary chancres by sex
P_VISIBLE_F = 0.30
P_VISIBLE_M = 0.50


def set_font(size=18):
    pl.rcParams.update({
        'font.size': size,
        'axes.titlesize': size,
        'axes.labelsize': size,
        'xtick.labelsize': size - 2,
        'ytick.labelsize': size - 2,
        'legend.fontsize': size - 4,
    })


def compute_stepdown(sw_prev_df, cs, start_year=2015, end_year=2025):
    """Compute 4-level step-down by age, sex, and SW group."""

    syph = sw_prev_df[sw_prev_df.disease == 'syph'].copy()
    skip = {'0-15', '65+', '50-65'}
    syph = syph[~syph.age.isin(skip)]

    # Aggregate seroprevalence from calib_stats (ever_exposed)
    mask = (cs.index >= start_year) & (cs.index <= end_year)
    ee_f = cs.loc[mask, ('epi_ts.syph_ever_exposed_f', '50%')].mean()
    ee_m = cs.loc[mask, ('epi_ts.syph_ever_exposed_m', '50%')].mean()
    ee_fsw = cs.loc[mask, ('epi_ts.syph_ever_exposed_fsw', '50%')].mean()

    # Compute sero/prevalence scaling ratios by sex/SW
    # Use 'Overall' sw_group to get the aggregate prevalence for scaling
    prev_f = syph[(syph.sex == 'Female') & (syph.sw_group == 'Overall')].prevalence.mean()
    prev_m = syph[(syph.sex == 'Male') & (syph.sw_group == 'Overall')].prevalence.mean()
    prev_fsw = syph[(syph.sex == 'Female') & (syph.sw_group == 'SW')].prevalence.mean()

    sero_ratio_f = max(ee_f / prev_f, 1.0) if prev_f > 0 else 1.0
    sero_ratio_m = max(ee_m / prev_m, 1.0) if prev_m > 0 else 1.0
    sero_ratio_fsw = max(ee_fsw / prev_fsw, 1.0) if prev_fsw > 0 else 1.0

    # Add computed columns
    results = []
    for _, row in syph.iterrows():
        is_female = row.sex == 'Female'
        is_sw = row.sw_group == 'SW'

        # Seroprevalence: scale total prevalence by sex/SW-specific ratio
        if is_female and is_sw:
            sero_ratio = sero_ratio_fsw
        elif is_female:
            sero_ratio = sero_ratio_f
        else:
            sero_ratio = sero_ratio_m

        sero = row.prevalence * sero_ratio
        active = row.active_prevalence
        primary = active * FRAC_PRIMARY
        p_vis = P_VISIBLE_F if is_female else P_VISIBLE_M
        visible = primary * p_vis

        results.append(dict(
            age=row.age,
            sex=row.sex,
            sw_group=row.sw_group,
            par_idx=row.par_idx,
            seroprevalence=sero * 100,
            active=active * 100,
            primary=primary * 100,
            visible=visible * 100,
        ))

    import pandas as pd
    return pd.DataFrame(results)


def plot_stepdown_panels(df):
    """
    2x2 panel figure showing the step-down across prevalence levels.
    Each panel has the same Fig 1A bar structure (3 groups × age groups)
    but shows a different prevalence level. Y-axes are independent so
    the small levels are readable, but the step-down is conveyed by the
    dramatic scale change across panels.
    """
    age_groups = ['15-20', '20-25', '25-30', '30-35', '35-50']
    groups = [
        ('Female', 'SW', 'Women in transactional sex', '#e377c2'),
        ('Female', 'Overall', 'All girls/women', '#ff9896'),
        ('Male', 'Overall', 'All boys/men', '#aec7e8'),
    ]

    levels = [
        ('seroprevalence', '(A) Seroprevalence (ever exposed)'),
        ('active', '(B) Active syphilis (primary + secondary)'),
        ('primary', '(C) Primary syphilis (ulcerative)'),
        ('visible', '(D) Visible chancre (presenting as GUD)'),
    ]

    dfm = df.groupby(['age', 'sex', 'sw_group']).mean(numeric_only=True).reset_index()

    fig, axes = pl.subplots(2, 2, figsize=(22, 14))
    axes = axes.flatten()

    n_groups = len(groups)
    bar_width = 0.25
    x_base = np.arange(len(age_groups))

    for pi, (level_key, title) in enumerate(levels):
        ax = axes[pi]

        for gi, (sex, sw_group, group_label, color) in enumerate(groups):
            x = x_base + (gi - 1) * bar_width
            heights = []
            for age in age_groups:
                sub = dfm[(dfm.age == age) & (dfm.sex == sex) & (dfm.sw_group == sw_group)]
                val = sub[level_key].values[0] if len(sub) > 0 else 0
                heights.append(val)
            ax.bar(x, heights, width=bar_width * 0.85, color=color, label=group_label if pi == 0 else None,
                   edgecolor='white', linewidth=0.5)

        ax.set_xticks(x_base)
        ax.set_xticklabels(age_groups)
        ax.set_title(title)
        ax.set_ylabel('Prevalence (%)')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        # Add y-axis max annotation to emphasize scale change
        ymax = ax.get_ylim()[1]
        ax.text(0.98, 0.95, f'max {ymax:.1f}%', transform=ax.transAxes,
                ha='right', va='top', fontsize=14, color='grey', style='italic')

    axes[0].legend(loc='upper left', framealpha=0.9)

    pl.tight_layout()
    return fig


if __name__ == '__main__':

    sw_prev_df = sc.loadobj(f'{RESULTS_DIR}/sw_prev_df.df')
    cs = sc.loadobj(f'{RESULTS_DIR}/zimbabwe_calib_stats_all.df')

    stepdown = compute_stepdown(sw_prev_df, cs)

    set_font(size=18)

    fig = plot_stepdown_panels(stepdown)

    pl.savefig(f'{FIGURES_DIR}/fig_stepdown.png', dpi=200, bbox_inches='tight')
    print(f'Saved {FIGURES_DIR}/fig_stepdown.png')
